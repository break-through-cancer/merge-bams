// #!/usr/bin/env nextflow
nextflow.enable.dsl = 2

import java.security.MessageDigest

/*
------------------------------------------------------------
Parameters + defaults
------------------------------------------------------------
*/
params.outdir      = params.outdir      ?: "results"
params.gatk_docker = params.gatk_docker ?: "broadinstitute/gatk:4.5.0.0"
params.max_sm_len  = params.max_sm_len  ?: 200


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
Mutect2 wants one sample identity in the merged BAM header.
We derive it from all BAM names.
------------------------------------------------------------
*/
def smParts = bamList.collect { p ->
    def n = new File(p.toString()).getName()
    n = n.replaceAll(/\.bam(\.bai)?$/, '')
    n = n.replaceAll(/\.sorted$/, '')
    n
}

def smJoined = smParts.join('__')

def unifiedSM =
    (smJoined.size() <= (params.max_sm_len as int))
        ? smJoined
        : "MERGED_${smParts.size()}_${md5(smJoined).substring(0,12)}"

println "Will merge ${bamList.size()} BAM(s)"
println "Unified SM = ${unifiedSM}"


/*
============================================================
PROCESS — MERGE_ALL_BAMS_AND_FIX_SM_HEADER
============================================================
What this does:
- Merge BAMs with MergeSamFiles
- Keep existing read groups / RGIDs
- Rewrite ONLY the SM tag in each @RG header line
- Reheader merged BAM
- Reindex final BAM
============================================================
*/
process MERGE_ALL_BAMS_AND_FIX_SM_HEADER {

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

    # provenance
    cat > merged.sm_manifest.txt << EOF
unified_SM=${unified_sm}
inputs:
${manifestText}
EOF

    echo "=== MERGING BAMs ==="
    gatk MergeSamFiles \\
        ${inputFlags} \\
        -O merged.preheader.bam

    echo "=== EXTRACT ORIGINAL HEADER ==="
    samtools view -H merged.preheader.bam > merged.header.sam

    echo "=== REWRITE ONLY SM ON @RG LINES ==="
    awk -v SM="${unified_sm}" 'BEGIN{FS=OFS="\\t"}
        /^@RG/ {
            found=0
            for (i=1; i<=NF; i++) {
                if (\$i ~ /^SM:/) {
                    \$i = "SM:" SM
                    found=1
                }
            }
            if (!found) {
                \$0 = \$0 OFS "SM:" SM
            }
            print
            next
        }
        { print }
    ' merged.header.sam > merged.header.smfixed.sam

    echo "=== VALIDATE HEADER SM VALUES ==="
    grep '^@RG' merged.header.smfixed.sam > rg_lines.txt || true

    if [[ ! -s rg_lines.txt ]]; then
        echo "ERROR: merged BAM header has no @RG lines; cannot define SM for Mutect2." >&2
        exit 1
    fi

    awk 'BEGIN{FS="\\t"}
        /^@RG/ {
            for (i=1; i<=NF; i++) {
                if (\$i ~ /^SM:/) {
                    sub(/^SM:/, "", \$i)
                    print \$i
                }
            }
        }
    ' merged.header.smfixed.sam | sort -u > sm_values.txt

    echo "SM values after rewrite:"
    cat sm_values.txt

    if [[ \$(wc -l < sm_values.txt) -ne 1 ]]; then
        echo "ERROR: expected exactly one SM after header rewrite" >&2
        exit 1
    fi

    if [[ "\$(cat sm_values.txt)" != "${unified_sm}" ]]; then
        echo "ERROR: final SM does not match expected unified SM" >&2
        exit 1
    fi

    echo "=== REHEADER MERGED BAM ==="
    samtools reheader -o merged.bam merged.header.smfixed.sam merged.preheader.bam

    echo "=== INDEX FINAL BAM ==="
    samtools index -o merged.bam.bai merged.bam

    test -f merged.bam
    test -f merged.bam.bai

    ls -lah merged.preheader.bam merged.bam merged.bam.bai merged.sm_manifest.txt
    """
}


/*
============================================================
WORKFLOW
============================================================
*/
workflow {

    Channel
        .fromList(bamList)
        .map { file(it) }
        .collect()
        .set { all_bams_ch }

    def manifest_lines = bamList.collect { "  - ${it}" }

    MERGE_ALL_BAMS_AND_FIX_SM_HEADER(
        all_bams_ch,
        unifiedSM,
        manifest_lines
    )
}