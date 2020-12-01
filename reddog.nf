#!/usr/bin/env nextflow
// Enable DSL2
nextflow.enable.dsl=2


// Import processes, workflows, channel helpers, and utility functions
// NOTE: imports separated for readability
// Alignment and variant calling processes
include { align_reads_se } from './src/processes/alignment.nf'
include { align_reads_pe } from './src/processes/alignment.nf'
include { call_snps } from './src/processes/variant_calling.nf'

// Mapping stats processes
include { calculate_gene_coverage_depth } from './src/processes/mapping_stats.nf'
include { calculate_mapping_statistics } from './src/processes/mapping_stats.nf'
include { aggregate_mapping_statistics } from './src/processes/mapping_stats.nf'
include { get_passing_replicons } from './src/processes/mapping_stats.nf'

// Allele matrix processes
include { create_allele_matrix } from './src/processes/allele_matrix.nf'
include { aggregate_allele_matrices } from './src/processes/allele_matrix.nf'
include { filter_allele_matrix } from './src/processes/allele_matrix.nf'

// Other processes
include { prepare_reference } from './src/processes/misc.nf'
include { subsample_reads_se } from './src/processes/misc.nf'
include { subsample_reads_pe } from './src/processes/misc.nf'
include { create_read_quality_report } from './src/processes/misc.nf'
include { create_mpileups } from './src/processes/misc.nf'
include { aggregate_snp_sites } from './src/processes/misc.nf'
include { aggregate_read_quality_reports } from './src/processes/misc.nf'
include { determine_coding_consequences } from './src/processes/misc.nf'
include { create_snp_alignment } from './src/processes/misc.nf'
include { infer_phylogeny } from './src/processes/misc.nf'

// Channel helper functions
include { get_read_prefix_and_type } from './src/channel_helpers.nf'
include { collect_passing_isolate_replicons } from './src/channel_helpers.nf'
include { sort_allele_matrices } from './src/channel_helpers.nf'
include { remove_empty_allele_matrices } from './src/channel_helpers.nf'

// Utility functions
include { print_splash } from './src/utilities.nf'
include { check_arguments } from './src/utilities.nf'
include { check_output_dir } from './src/utilities.nf'
include { check_disallowed_arguments } from './src/utilities.nf'
include { check_host } from './src/utilities.nf'
include { check_boolean_option } from './src/utilities.nf'
include { write_reference_data_to_run_config } from './src/utilities.nf'
include { write_param_data_to_run_config } from './src/utilities.nf'
include { validate_merge_data } from './src/utilities.nf'

// Workflows
include { merge } from './src/merge_workflow.nf'


// Check configuration
print_splash()
check_disallowed_arguments(workflow)
check_arguments(params)
check_output_dir(params)
check_host(workflow)


// Require some variables to be boolean
// We must check and change values if needed
run_read_subsample = check_boolean_option(params.subsample_reads, 'subsample_reads')
run_read_quality_report = check_boolean_option(params.read_quality_report, 'read_quality_report')
run_phylogeny = check_boolean_option(params.force_tree, 'force_tree')
run_merge = check_boolean_option(params.merge_run, 'merge_run')
merge_ignore_errors = check_boolean_option(params.merge_ignore_errors, 'merge_ignore_errors')


// Check integer params
if (run_read_subsample & params.subsample_read_count.getClass() != java.lang.Integer) {
  if (! params.subsample_read_count.isInteger()) {
    exit 1, "ERROR: subsample_read_count must be an integer, got {params.subsample_read_count}"
  }
}


// Get reads and seperate into pe and se channels based on prefix
reads = Channel.fromPath(params.reads).ifEmpty {
    exit 1, "ERROR: did not find any read files with '${params.reads}'"
  }.map {
    get_read_prefix_and_type(it)
  }.branch {
    paired: it[0] == 'pe'
    single: it[0] == 'se'
}
reads_se = reads.single.map { it[1..-1] }.groupTuple()
reads_pe = reads.paired.map { it[1..-1] }.groupTuple()
// Check that we have the expected number of reads for each prefix in pe and se channels and flatten tuple
reads_pe = reads_pe.map {
  if (it[1].size() != 2) {
    exit 1, "ERROR: didn't get exactly two readsets prefixed with ${it[0]}:\n${it[1]}"
  }
  [it[0], *it[1]]
}
reads_se = reads_se.map {
  if (it[1].size() != 1) {
    exit 1, "ERROR: didn't get exactly one readset prefixed with ${it[0]}:\n${it[1]}"
  }
  [it[0], *it[1]]
}


// TODO: check we have correct number of reads in output


// Create file object for reference and check it exists
reference_fp = file(params.reference)
if (! reference_fp.exists()) {
  exit 1, "ERROR: reference input '${reference_fp}' does not exist"
}


// Create run config output file
write_reference_data_to_run_config()
write_param_data_to_run_config()


// Prepare for a merge run if requested
if (run_merge) {
  // Validate data to be merged
  if (params.existing_run_dir.isEmpty()) {
    exit 1, "ERROR: a merge run requires existing_run_dir to be set"
  }
  merge_source_dir = file(params.existing_run_dir)
  if (! merge_source_dir.exists()) {
    exit 1, "ERROR: directory for existing_run_dir (${params.existing_run_dir}) does not exist"
  }
  validate_merge_data(merge_ignore_errors)
  // Create channels for merge data
  merge_source_bams = Channel.fromPath(merge_source_dir / 'bams/*.bam')
  merge_source_vcfs = Channel.fromPath(merge_source_dir / 'vcfs/*.vcf')
  merge_source_fastqc = Channel.fromPath(merge_source_dir / 'fastqc/individual_reports/*{zip,html}')
  merge_source_gene_depth = Channel.fromPath(merge_source_dir / '*gene_depth.csv')
  merge_source_gene_coverage = Channel.fromPath(merge_source_dir / '*gene_coverage.csv')
  merge_source_mapping_stats = Channel.fromPath(merge_source_dir / '*mapping_stats.tsv')
  merge_source_allele_matrices = Channel.fromPath(merge_source_dir / '*alleles.csv')
}


log.info('----------------------------------------------------------------------------------')


workflow {
  main:
    reference_data = prepare_reference(reference_fp)

    if (run_read_subsample) {
      reads_se = subsample_reads_se(reads_se)
      reads_pe = subsample_reads_pe(reads_pe)
    }

    ch_fastqc = Channel.empty()
    if (run_read_quality_report) {
      // Create flat channel containing only readset filepaths
      reads_all_fps = reads_pe.map { it[1..-1] }.mix(reads_se.map { it[1..-1] }).flatten()
      ch_fastqc = create_read_quality_report(reads_all_fps).output
    }

    align_data_se = align_reads_se(reads_se, reference_data.fasta, reference_data.bt2_index)
    align_data_pe = align_reads_pe(reads_pe, reference_data.fasta, reference_data.bt2_index)
    align_data_bams = align_data_se.bams.mix(align_data_pe.bams)
    align_data_bams_unmapped = align_data_se.bams_unmapped.mix(align_data_pe.bams_unmapped)

    mpileup_data = create_mpileups(align_data_bams)

    bams_and_mpileups = align_data_bams.join(mpileup_data.output)
    snp_data = call_snps(bams_and_mpileups, reference_data.fasta, reference_data.samtools_index)

    bams_vcfs_and_stats = align_data_bams.join(align_data_bams_unmapped)
      .join(snp_data.vcfs)
      .join(snp_data.coverage_depth)

    gene_coverage_depth = calculate_gene_coverage_depth(mpileup_data.output_noid.collect(), reference_fp)

    stats_data = calculate_mapping_statistics(bams_vcfs_and_stats)

    stats_aggregated_data = aggregate_mapping_statistics(stats_data.collect(), reference_data.name)

    // These channels are created so that if we're doing a merge run they will be modified inplace
    ch_mapping_stats = stats_aggregated_data.stats
    ch_snp_sites = snp_data.sites
    ch_bams = align_data_bams

    // Merge existing run data if requested
    if (run_merge) {
      // NOTE: reference_fp.simpleName is used over reference_data.name as we need it immediately in some functions
      merge_data = merge(
        merge_source_bams,
        merge_source_vcfs,
        ch_fastqc,
        merge_source_fastqc,
        gene_coverage_depth.depth,
        merge_source_gene_depth,
        gene_coverage_depth.coverage,
        merge_source_gene_coverage,
        stats_aggregated_data.stats,
        merge_source_mapping_stats,
        merge_source_allele_matrices,
        run_read_quality_report,
        reference_fp.simpleName,
      )
      ch_fastqc = merge_data.fastqc
      ch_snp_sites = ch_snp_sites.mix(merge_data.snp_sites)
      ch_mapping_stats = merge_data.mapping_stats
      ch_bams = ch_bams.mix(merge_data.bams)
    }

    if (run_read_quality_report) {
      aggregate_read_quality_reports(ch_fastqc.collect())
    }

    isolate_replicon_data = get_passing_replicons(ch_mapping_stats.collect())

    sites_data = aggregate_snp_sites(ch_snp_sites.collect(), isolate_replicon_data.output)

    // For each isolate get replicons that passing mapping stats criteria and have at least one SNP
    // Filter isolates that do not pass
    isolate_replicons_passing = collect_passing_isolate_replicons(isolate_replicon_data.output)

    bams_and_sites = ch_bams.join(isolate_replicons_passing).combine(sites_data.output)
    allele_matrix_data = create_allele_matrix(bams_and_sites, reference_data.fasta)

    allele_matrices_by_replicon = sort_allele_matrices(allele_matrix_data.output)
    allele_matrix_aggregate_data = aggregate_allele_matrices(allele_matrices_by_replicon, sites_data.output, reference_data.name)
    allele_matrix_aggregate_data = remove_empty_allele_matrices(allele_matrix_aggregate_data.output)

    allele_matrices_cons = filter_allele_matrix(allele_matrix_aggregate_data, reference_data.name)

    determine_coding_consequences(allele_matrices_cons, reference_fp)

    snp_alignment = create_snp_alignment(allele_matrices_cons, reference_data.name)
    infer_phylogeny(snp_alignment, reference_data.name, isolate_replicons_passing.count(), run_phylogeny)
}
