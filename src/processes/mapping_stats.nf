// Get gene coverage and mean depth
process calculate_gene_coverage_depth {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }, mode: 'copy'

  input:
  path mpileup_fps
  path reference_gbk_fp

  output:
  path 'gene_coverage.csv', emit: coverage
  path 'gene_depth.csv', emit: depth

  script:
  reference_name = reference_gbk_fp.simpleName
  """
  create_coverage_depth_matrices.py --mpileup_fps ${mpileup_fps} --reference_fp ${reference_gbk_fp} --output_dir ./
  """
}


// Calculate mapping statistics
//   - get unmapped reads from unmapped reads BAM
//   - coverage and mapped read counts previously calculated from mpileups
//   - get SNP, INDEL, and heterzygous SNP counts from VCFs
//   - calculate appropriate statisitics
//   - fail isolates on depth or coverage
//   - additionally for the largest replicon, require that at least some proportion of reads mapped
process calculate_mapping_statistics {
  input:
  tuple val(isolate_id), path(bam_fp), path(bam_index_fp), path(bam_unmapped_fp), path(vcf_q30_fp), path(vcf_hets_fp), path(coverage_depth_fp)

  output:
  path "${isolate_id}_mapping_stats.tsv", emit: output

  script:
  """
  calculate_mapping_stats.py --bam_fp ${bam_fp} --vcf_q30_fp ${vcf_q30_fp} --vcf_hets_fp ${vcf_hets_fp} \
    --coverage_depth_fp ${coverage_depth_fp} --bam_unmapped_fp ${bam_unmapped_fp} \
    --min_depth ${params.mapping_depth_min} --min_coverage ${params.mapping_cover_min} \
    --min_mapped_reads ${params.mapping_mapped_min} > ${isolate_id}_mapping_stats.tsv
  """
}


// Aggregate mapping statistics
//   - read in all records and aggregate by replicon
//   - for passing isolates, calculate phylogeny group using the ratio of 'total reads mapped' : 'coverage' / 100
//     - set maximum ratio for ingroup as 'ratio mean' + 'ratio stddev' * 'modifier specified by user'
//   - write aggregate stats to file, in the order of 'fail', 'outgroup', 'ingroup'
//   - get a list of isolates that pass mapping on at least one replicon
//   - get a list of replicons that have any passing isolate
//     - determined by having pass status and more than one SNP
process aggregate_mapping_statistics {
  publishDir "${params.output_dir}", pattern: '*_mapping_stats.tsv', saveAs: { filename -> "${reference_name}_${filename}" }, mode: 'copy'

  input:
  path(mapping_stats_fps)
  val reference_name

  output:
  path '*_mapping_stats.tsv', emit: stats

  script:
  """
  aggregate_mapping_stats.py --rep_stats_fps ${mapping_stats_fps} --output_dir ./ --stddev_mod ${params.outgroup_mod}
  """
}


process get_passing_replicons {
  executor 'local'

  input:
  path stats_fps

  output:
  path 'passing_isolate_replicons.tsv', emit: output

  script:
  """
  get_passing_isolate_replicons.awk ${stats_fps} > passing_isolate_replicons.tsv
  """
}
