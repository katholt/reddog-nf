process merge_gene_depth {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  path file1, stageAs: 'file_1.tsv'
  path file2, stageAs: 'file_2.tsv'

  output:
  path 'gene_depth.tsv'

  script:
  """
  touch gene_depth.tsv
  #./merge_gene_stat_matrix.py --fp_1 file_1.tsv --fp_2 file_2.tsv > gene_depth.tsv
  """
}


process merge_gene_coverage {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  path file1, stageAs: 'file_1.tsv'
  path file2, stageAs: 'file_2.tsv'

  output:
  path 'gene_coverage.tsv'

  script:
  """
  touch gene_coverage.tsv
  #./merge_gene_stat_matrix.py --fp_1 file_1.tsv --fp_2 file_2.tsv > gene_coverage.tsv
  """
}


process merge_allele_matrices {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${replicon_id}_alleles.tsv" }

  input:
  tuple replicon_id, path(files)

  output:
  path 'output.tsv'

  script:
  file1 = files[0]
  file2 = files[1]
  """
  touch output.tsv
  #./merge_allele_matrices.py --fp_1 ${file1} --fp_2 ${file2} > output.tsv
  """
}


process merge_mapping_stats {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${replicon_id}_mapping_stats.tsv" }

  input:
  tuple replicon_id, path(files)

  output:
  path 'output.tsv'

  script:
  file1 = files[0]
  file2 = files[1]
  """
  touch output.tsv
  #./merge_mapping_stats.py --fp_1 ${file1} --fp_2 ${file2} > output.tsv
  """
}
