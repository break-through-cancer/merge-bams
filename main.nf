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
  """
  set -euo pipefail

  echo "=== MERGE_ALL_BAMS ==="
  echo "Inputs:"
  ls -lah

  gatk MergeSamFiles \\
    $(for b in ${bams}; do echo -n " -I \$b"; done) \\
    -O merged.bam \\
    --CREATE_INDEX true

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