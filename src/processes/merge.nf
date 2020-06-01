process gene_depth {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  path file1, stageAs: 'file_1.tsv'
  path file2, stageAs: 'file_2.tsv'
  val reference_name

  output:
  path 'gene_depth.tsv'

  script:
  """
  merge_gene_stat_matrix.py --fp_1 file_1.tsv --fp_2 file_2.tsv > gene_depth.tsv
  """
}


process gene_coverage {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  path file1, stageAs: 'file_1.tsv'
  path file2, stageAs: 'file_2.tsv'
  val reference_name

  output:
  path 'gene_coverage.tsv'

  script:
  """
  merge_gene_stat_matrix.py --fp_1 file_1.tsv --fp_2 file_2.tsv > gene_coverage.tsv
  """
}


process allele_matrix {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${replicon_id}_alleles.tsv" }

  input:
  tuple replicon_id, path(files)
  val reference_name

  output:
  tuple val(replicon_id), path('output.tsv')

  script:
  file1 = files[0]
  file2 = files[1]
  """
  merge_allele_matrices.py --fp_1 ${file1} --fp_2 ${file2} > output.tsv
  """
}


process mapping_stats {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${replicon_id}_mapping_stats.tsv" }

  input:
  tuple replicon_id, path(files)
  val reference_name

  output:
  path 'output.tsv'

  script:
  file1 = files[0]
  file2 = files[1]
  """
  merge_mapping_stats.py --fp_1 ${file1} --fp_2 ${file2} --stddev_mod ${params.outgroup_mod} > output.tsv
  """
}
