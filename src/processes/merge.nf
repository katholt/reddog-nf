// Merge and publish fastqc report dataset
process merge_fastqc_data {
  publishDir "${params.output_dir}/fastqc/individual_reports/"

  input:
  path fastqc_fp

  output:
  path '*html'
  path '*zip', emit: output

  script:
  """
  fastqc -o . -t 1 -q -f fastq ${reads_fwd} ${reads_rev}
  """
}

