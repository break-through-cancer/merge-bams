#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

import java.security.MessageDigest

// #!/usr/bin/env nextflow
// nextflow.enable.dsl = 2

// params.outdir      = params.outdir      ?: "results"
// params.gatk_docker = params.gatk_docker ?: "broadinstitute/gatk:4.5.0.0"


// def bamList = params.bams

// if( !bamList || bamList.isEmpty() ) {
//   error "Missing required --bams"
// }

// process MERGE_ALL_BAMS {
//   tag "merge_all"
//   container "${params.gatk_docker}"
//   publishDir "${params.outdir}", mode: 'copy'

//   input:
//     path bams

//   output:
//     path "merged.bam"
//     path "merged.bam.bai"

//   script:
//     def inputFlags = bams.collect { "-I ${it}" }.join(' ')
//     """
//     set -euo pipefail

//     gatk MergeSamFiles \
//       ${inputFlags} \
//       -O merged.bam \
//       --CREATE_INDEX true

//     # GATK/Picard may output merged.bai (not merged.bam.bai)
//     if [[ -f merged.bai && ! -f merged.bam.bai ]]; then
//       mv merged.bai merged.bam.bai
//     fi

//     # If still missing, try to create one (if samtools exists in image)
//     if [[ ! -f merged.bam.bai ]]; then
//       if command -v samtools >/dev/null 2>&1; then
//         samtools index -o merged.bam.bai merged.bam
//       fi
//     fi

//     # Hard fail if index truly not produced
//     test -f merged.bam.bai

//     ls -lah merged.bam merged.bam.bai
//     """
// }


// workflow {
//   Channel
//     .fromList(bamList)
//     .map { file(it) }
//     .collect()
//     .set { all_bams_ch }

//   MERGE_ALL_BAMS(all_bams_ch)
// }


/*
Read-group defaults
RGPL = sequencing platform
RGLB = library name
*/
params.rgpl        = params.rgpl ?: "ILLUMINA"
params.rglb        = params.rglb ?: "lib1"

/*
Maximum allowed SM length before hashing
(Some tools break with extremely long SMs)
*/
params.max_sm_len  = params.max_sm_len ?: 200


/*
------------------------------------------------------------
Input BAM handling
------------------------------------------------------------
*/
def bamList = params.bams

if( !bamList || (bamList instanceof List && bamList.isEmpty()) )
    error "Missing required --bams"

/* allow comma-separated CLI input */
if( !(bamList instanceof List) )
    bamList = bamList.toString().split(/\s*,\s*/).findAll{ it }


/*
------------------------------------------------------------
Helper: MD5 hash generator
Used only if merged sample name becomes too long
------------------------------------------------------------
*/
String md5(String s) {
    MessageDigest
        .getInstance("MD5")
        .digest(s.getBytes("UTF-8"))
        .encodeHex()
        .toString()
}


/*
------------------------------------------------------------
Create ONE unified sample name (SM)
------------------------------------------------------------
Mutect2 requires exactly ONE SM.
We derive it from all BAM names.
*/
def smParts = bamList.collect { p ->
    def n = new File(p.toString()).getName()
    n = n.replaceAll(/\.bam(\.bai)?$/, '')
    n = n.replaceAll(/\.sorted$/, '')
    n
}

def smJoined = smParts.join('__')

def unifiedSM =
    (smJoined.size() <= params.max_sm_len)
        ? smJoined
        : "MERGED_${smParts.size()}_${md5(smJoined).substring(0,12)}"

println "Will merge ${bamList.size()} BAM(s)"
println "Unified SM = ${unifiedSM}"


/*
============================================================
PROCESS 1 — FIX_READGROUPS
============================================================

Why this exists:

 ONE SM,  UNIQUE RGID,  RGPU
*/
process FIX_READGROUPS {

    tag { "fix_rg_${idx}" }
    container "${params.gatk_docker}"

    input:
        tuple path(bam), val(idx)
        val unified_sm
        val rgpl
        val rglb

    output:
        tuple path("fixed.${idx}.bam"),
              path("fixed.${idx}.bam.bai"),
              val(idx),
              emit: fixed

    script:
    """
    set -euo pipefail

    # Unique identifiers per BAM
    RGID="rg${idx}"
    RGPU="pu${idx}"

    gatk AddOrReplaceReadGroups \\
        -I "${bam}" \\
        -O "fixed.${idx}.bam" \\
        -RGID "\$RGID" \\
        -RGLB "${rglb}" \\
        -RGPL "${rgpl}" \\
        -RGPU "\$RGPU" \\
        -RGSM "${unified_sm}" \\
        --CREATE_INDEX true

    # normalize index naming
    if [[ -f fixed.${idx}.bai && ! -f fixed.${idx}.bam.bai ]]; then
        mv fixed.${idx}.bai fixed.${idx}.bam.bai
    fi

    # guarantee index exists
    if [[ ! -f fixed.${idx}.bam.bai ]]; then
        gatk BuildBamIndex \
            -I fixed.${idx}.bam \
            -O fixed.${idx}.bam.bai
    fi
    """
}



/*
============================================================
PROCESS 2 — MERGE_ALL_BAMS
============================================================
Merge already-normalized BAMs safely
*/
process MERGE_ALL_BAMS {

    tag "merge_all"
    container "${params.gatk_docker}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path bams
        val unified_sm
        val manifest_lines

    output:
        path "merged.bam"
        path "merged.bam.bai"
        path "merged.sm_manifest.txt"

    script:
        def inputFlags = bams.collect { "-I ${it}" }.join(' ')
        def manifestText = manifest_lines.join('\n')

    """
    set -euo pipefail

    # record provenance
    cat > merged.sm_manifest.txt << EOF
unified_SM=${unified_sm}
inputs:
${manifestText}
EOF

    gatk MergeSamFiles \\
        ${inputFlags} \\
        -O merged.bam \\
        --CREATE_INDEX true

    # normalize index naming
    if [[ -f merged.bai && ! -f merged.bam.bai ]]; then
        mv merged.bai merged.bam.bai
    fi

    if [[ ! -f merged.bam.bai ]]; then
        gatk BuildBamIndex \
            -I merged.bam \
            -O merged.bam.bai
    fi

    test -f merged.bam.bai
    """
}



/*
============================================================
WORKFLOW
============================================================
*/
workflow {

    /*
    enumerate() replaces withIndex()
    (Cirro-safe)
    */
    Channel
        .fromList(bamList)
        .map { file(it) }
        .enumerate()
        .map { idx, bam -> tuple(bam, idx) }
        .set { bam_idx_ch }


    /*
    Step 1: normalize read groups
    */
    FIX_READGROUPS(
        bam_idx_ch,
        unifiedSM,
        params.rgpl,
        params.rglb
    )


    /*
    Step 2: collect fixed BAMs
    */
    FIX_READGROUPS.out.fixed
        .map { fixed_bam, fixed_bai, idx -> fixed_bam }
        .collect()
        .set { fixed_bams_ch }


    /*
    Step 3: merge
    */
    def manifest_lines =
        bamList.collect { "  - ${it}" }

    MERGE_ALL_BAMS(
        fixed_bams_ch,
        unifiedSM,
        manifest_lines
    )
}