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

#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

import java.security.MessageDigest

params.outdir      = params.outdir      ?: "results"
params.gatk_docker = params.gatk_docker ?: "broadinstitute/gatk:4.5.0.0"
params.rgpl        = params.rgpl        ?: "ILLUMINA"
params.rglb        = params.rglb        ?: "lib1"
params.max_sm_len  = params.max_sm_len  ?: 200   // keep SM reasonably sized

def bamList = params.bams
if( !bamList || (bamList instanceof List && bamList.isEmpty()) ) {
  error "Missing required --bams"
}
if( !(bamList instanceof List) ) {
  // allow comma-separated string
  bamList = bamList.toString().split(/\s*,\s*/).findAll{ it }
}

String md5(String s) {
  MessageDigest.getInstance("MD5").digest(s.getBytes("UTF-8")).encodeHex().toString()
}

// Build a unified SM that reflects all BAMs.
// If it gets too long, fall back to a hashed name but still deterministic.
def smParts = bamList.collect { p ->
  def n = new File(p.toString()).getName()
  n = n.replaceAll(/\.bam(\.bai)?$/, '')
  n = n.replaceAll(/\.sorted$/, '')
  return n
}
def smJoined = smParts.join('__')
def unifiedSM = (smJoined.size() <= (params.max_sm_len as int))
  ? smJoined
  : "MERGED_${smParts.size()}_${md5(smJoined).substring(0,12)}"

println "Will merge ${bamList.size()} BAM(s)"
println "Unified SM: ${unifiedSM}"

process FIX_READGROUPS {
  tag { "fix_rg_${idx}" }
  container "${params.gatk_docker}"

  input:
    tuple path(bam), val(idx)
    val unified_sm
    val rgpl
    val rglb

  output:
    tuple path("fixed.${idx}.bam"), path("fixed.${idx}.bam.bai"), val(idx), emit: fixed

  script:
  """
  set -euo pipefail

  # Unique per-input RG identifiers (avoid collisions when merging)
  RGID="rg${idx}"
  RGPU="pu${idx}"
  RGLB="${rglb}"
  RGPL="${rgpl}"
  RGSM="${unified_sm}"

  gatk AddOrReplaceReadGroups \\
    -I "${bam}" \\
    -O "fixed.${idx}.bam" \\
    -RGID "\${RGID}" \\
    -RGLB "\${RGLB}" \\
    -RGPL "\${RGPL}" \\
    -RGPU "\${RGPU}" \\
    -RGSM "\${RGSM}" \\
    --CREATE_INDEX true || true

  # Some versions write fixed.\${idx}.bai instead of fixed.\${idx}.bam.bai; normalize
  if [[ -f "fixed.${idx}.bai" && ! -f "fixed.${idx}.bam.bai" ]]; then
    mv "fixed.${idx}.bai" "fixed.${idx}.bam.bai"
  fi

  # If still missing, build an index
  if [[ ! -f "fixed.${idx}.bam.bai" ]]; then
    gatk BuildBamIndex -I "fixed.${idx}.bam" -O "fixed.${idx}.bam.bai" || true
  fi

  test -f "fixed.${idx}.bam.bai"
  """
}

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

    # Record what SM we enforced and what we merged
    cat > merged.sm_manifest.txt << 'EOF'
unified_SM=${unified_sm}
inputs:
${manifestText}
EOF

    gatk MergeSamFiles \\
      ${inputFlags} \\
      -O merged.bam \\
      --CREATE_INDEX true

    # Normalize index name (merged.bai vs merged.bam.bai)
    if [[ -f merged.bai && ! -f merged.bam.bai ]]; then
      mv merged.bai merged.bam.bai
    fi

    if [[ ! -f merged.bam.bai ]]; then
      gatk BuildBamIndex -I merged.bam -O merged.bam.bai || true
    fi

    test -f merged.bam.bai
    ls -lah merged.bam merged.bam.bai merged.sm_manifest.txt
    """
}

workflow {
  // Make tuples (bam, idx)
  Channel
    .fromList(bamList)
    .map { file(it) }
    .withIndex()
    .map { bam, idx -> tuple(bam, idx) }
    .set { bam_idx_ch }

  // Rewrite RGs to unify SM + ensure unique RGID/RGPU
  FIX_READGROUPS(
    bam_idx_ch,
    unifiedSM,
    params.rgpl as String,
    params.rglb as String
  )

  // Collect fixed BAMs for merge
  def fixed_bams_ch = FIX_READGROUPS.out.fixed
    .map { fixed_bam, fixed_bai, idx -> fixed_bam }
    .collect()

  // Build manifest lines (original names) for output sidecar
  def manifest_lines = bamList.collect { "  - ${it}" }

  MERGE_ALL_BAMS(
    fixed_bams_ch,
    unifiedSM,
    manifest_lines
  )
}
