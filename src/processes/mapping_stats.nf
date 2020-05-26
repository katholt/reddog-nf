// Get gene coverage and mean depth
process calculate_gene_coverage_depth {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  path mpileup_fps
  path reference_gbk_fp

  output:
  path 'gene_coverage.tsv', emit: coverage
  path 'gene_depth.tsv', emit: depth

  script:
  reference_name = reference_gbk_fp.simpleName
  """
  create_coverage_depth_matrices.py --mpileup_fps ${mpileup_fps} --reference_fp ${reference_gbk_fp} --output_dir ./
  """
}


// Calculate mapping statistics
//   - get total reads from mapping metrics file
//   - coverage and mapped read counts previously calculated from mpileups
//   - get SNP, INDEL, and heterzygous SNP counts from VCFs
//   - calculate appropriate statisitics
//   - fail isolates on depth or coverage
//   - additionally for the largest replicon, require that at least 50% reads mapped
process calculate_mapping_statistics {
  input:
  tuple isolate_id, path(bam_fp), path(bam_index_fp), path(vcf_q30_fp), path(vcf_hets_fp), path(coverage_depth_fp), path(mapping_metrics_fp)

  output:
  path "${isolate_id}_mapping_stats.tsv", emit: output

  script:
  """
  calculate_mapping_stats.py --bam_fp ${bam_fp} --vcf_q30_fp ${vcf_q30_fp} --vcf_hets_fp ${vcf_hets_fp} --coverage_depth_fp ${coverage_depth_fp} --mapping_metrics_fp ${mapping_metrics_fp} > ${isolate_id}_mapping_stats.tsv
  """
}


// Aggregate mapping statistics
//   - read in all records and aggregate by replicon
//   - for passing isolates, calculate phylogeny group using the ratio of 'total reads mapped' : 'coverage' / 100
//     - set maximum ratio for ingroup as 'ratio mean' + 'ratio stddev' * 2
//   - write aggregate stats to file, in the order of 'fail', 'outgroup', 'ingroup'
//   - get a list of isolates that pass mapping on at least one replicon
//   - get a list of replicons that have any passing isolate
//     - determined by having pass status and more than one SNP
process aggregate_mapping_statistics {
  publishDir "${params.output_dir}", pattern: '*_mapping_stats.tsv', saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  path(mapping_stats_fps)
  val reference_name

  output:
  path '*_mapping_stats.tsv', emit: stats
  path 'isolate_replicons_passing.tsv', emit: isolate_replicons

  script:
  """
  aggregate_mapping_stats.py --rep_stats_fps ${mapping_stats_fps} --output_dir ./
  get_passing_isolate_replicons.awk ${mapping_stats_fps} > isolate_replicons_passing.tsv
  """
}



