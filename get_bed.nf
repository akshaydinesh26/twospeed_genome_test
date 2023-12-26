#!/usr/local/bin/nextflow

params.gff_file="$projectDir/gff/*.gff"
parmas.outdir="$projectDir/bedfiles"

process processGff {

  publishDir "${params.outdir}", mode: 'copy'
  input:
  path gff_file
  tuple val($sample_id), path(gff_file)
  
  output:
  path "$sample_id"

  script:
  """
  ./prepare_gff.sh $gff_file $sample_id
  """
}

channel.fromPath(params.gff_file, checkIfExists: true).set{ input_ch }
workflow {
  processgff_ch = processGff(input_ch)
}