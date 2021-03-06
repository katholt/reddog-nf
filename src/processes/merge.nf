process copy_bams {
  publishDir "${params.output_dir}/bams/", mode: 'copy'

  executor 'local'

  input:
  path src_fn

  output:
  path src_fn

  script:
  """
  """
}

process copy_vcfs {
  publishDir "${params.output_dir}/vcfs/", mode: 'copy'

  executor 'local'

  input:
  path src_fn

  output:
  path src_fn

  script:
  """
  """
}

process copy_fastqc {
  publishDir "${params.output_dir}/fastqc/individual_reports/", mode: 'copy'

  executor 'local'

  input:
  path src_fn

  output:
  path src_fn

  script:
  """
  """
}

process index_bam {
  input:
  path bam_fp

  output:
  tuple val(isolate_id), path(bam_fp), path("*.bam.*"), emit: output

  script:
  isolate_id = bam_fp.simpleName
  """
  samtools index ${bam_fp}
  """
}

process gene_depth {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }, mode: 'copy'

  input:
  path file1, stageAs: 'file_1.csv'
  path file2, stageAs: 'file_2.csv'
  val reference_name

  output:
  path 'gene_depth.csv'

  script:
  """
  merge_gene_stat_matrix.py --fp_1 file_1.csv --fp_2 file_2.csv > gene_depth.csv
  """
}


process gene_coverage {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }, mode: 'copy'

  input:
  path file1, stageAs: 'file_1.csv'
  path file2, stageAs: 'file_2.csv'
  val reference_name

  output:
  path 'gene_coverage.csv'

  script:
  """
  merge_gene_stat_matrix.py --fp_1 file_1.csv --fp_2 file_2.csv > gene_coverage.csv
  """
}


process mapping_stats {
  executor 'local'

  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${replicon_id}_mapping_stats.tsv" }, mode: 'copy'

  input:
  tuple val(replicon_id), path(files)
  val reference_name

  output:
  path "${replicon_id}_mapping_stats.tsv", emit: output

  script:
  file1 = files[0]
  file2 = files[1]
  """
  # Rename symlinks so we can write to desired filename
  mv ${file1} file1.tsv
  mv ${file2} file2.tsv
  merge_mapping_stats.py --fp_1 file1.tsv --fp_2 file2.tsv --stddev_mod ${params.outgroup_mod} > ${replicon_id}_mapping_stats.tsv
  """
}


process collect_snp_sites {
  executor 'local'

  input:
  tuple val(replicon_id), path(allele_matrix)

  output:
  path "merge_${replicon_id}_sites.tsv"

  script:
  """
  echo -e "#CHROM\tPOS" > merge_${replicon_id}_sites.tsv
  tail -n+2 ${allele_matrix} | cut -f1 -d, | awk '{print "${replicon_id}\t"\$0 }' >> merge_${replicon_id}_sites.tsv
  """
}
