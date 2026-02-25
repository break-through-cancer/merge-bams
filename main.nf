#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.outdir      = params.outdir      ?: "results"
params.gatk_docker = params.gatk_docker ?: "broadinstitute/gatk:4.5.0.0"


def bamList = params.bams

if( !bamList || bamList.isEmpty() ) {
  error "Missing required --bams"
}

process MERGE_ALL_BAMS {
  tag "merge_all"
  container "${params.gatk_docker}"
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path bams

  output:
    path "merged.bam"
    path "merged.bam.bai"

  script:
    def inputFlags = bams.collect { "-I ${it}" }.join(' ')
    """
    set -euo pipefail

    gatk MergeSamFiles \
      ${inputFlags} \
      -O merged.bam \
      --CREATE_INDEX true

    # GATK/Picard may output merged.bai (not merged.bam.bai)
    if [[ -f merged.bai && ! -f merged.bam.bai ]]; then
      mv merged.bai merged.bam.bai
    fi

    # If still missing, try to create one (if samtools exists in image)
    if [[ ! -f merged.bam.bai ]]; then
      if command -v samtools >/dev/null 2>&1; then
        samtools index -o merged.bam.bai merged.bam
      fi
    fi

    # Hard fail if index truly not produced
    test -f merged.bam.bai

    ls -lah merged.bam merged.bam.bai
    """
}


workflow {
  Channel
    .fromList(bamList)
    .map { file(it) }
    .collect()
    .set { all_bams_ch }

  MERGE_ALL_BAMS(all_bams_ch)
}