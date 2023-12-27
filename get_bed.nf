#!/usr/local/bin/nextflow

params.gff_file="$projectDir/gff/*.gff"
params.outdir="$projectDir/bedfiles"
params.finaldir="$projectDir/finalFIR"

process processGff {

  publishDir "${params.outdir}", mode: 'copy'
  input:
  path gff_file

  output:
  path "${gff_file.baseName}"

  script:
  """
  echo $gff_file
  $projectDir/prepare_gff.sh $gff_file "${gff_file.baseName}"
  """
}

process firmaker {

   publishDir "${params.finaldir}", mode: 'copy'

  input:
  path bedfile

  output:
  path "${bedfile.baseName}"

  script:
  """
  Rscript "$projectDir/get_flanking_fir.R" $bedfile "${bedfile.baseName}" 
    
  """
}

channel.fromPath(params.gff_file, checkIfExists: true).set{ input_ch }
workflow {
   processGff_ch = processGff(input_ch)
   obtain_fir = firmaker(processGff_ch)
   obtain_fir.view()
}
