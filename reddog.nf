#!/usr/bin/env nextflow
// Enable DSL2
nextflow.preview.dsl=2


// Workflows
include preprocess_reads from './src/modules/read_preprocessing.nf'
include call_variants from './src/modules/variant_calling.nf'
include run_mapping_stats from './src/modules/mapping_stats.nf'
include run_allele_matrices from './src/modules/allele_matrices.nf'
include run_merge_results from './src/modules/merge_results.nf'
include run_post from './src/modules/post_analysis.nf'

// Processes used in main workflow - separated for readability
include prepare_reference from './src/processes/common.nf'
include { create_snp_alignment; infer_phylogeny } from './src/processes/common.nf'

// Utility - separated for readability
include { print_splash; print_help } from './src/utilities.nf'
include { check_arguments; check_input_files; check_output_dir } from './src/utilities.nf'
include { check_host; check_boolean_option } from './src/utilities.nf'


// Check configuration
print_splash()
check_arguments(params)
check_input_files(workflow, params)
check_output_dir(params)
check_host(workflow)


// Require optional stage variables to be boolean
// We must check and change values if needed. The global param variables are immutable so instead we declare new ones
run_read_subsample = check_boolean_option(params.subsample_reads, 'subsample_reads')
run_quality_assessment = check_boolean_option(params.quality_assessment, 'quality_assessment')
run_phylogeny = check_boolean_option(params.force_tree, 'force_tree')
merge_run = check_boolean_option(params.merge_run, 'merge_run')


// Check integer params
if (run_read_subsample & params.subsample_read_count.getClass() != java.lang.Integer) {
  if (! params.subsample_read_count.isInteger()) {
    exit 1, "error: subsample_read_count must be an integer, got {params.subsample_read_count}"
  }
}


// Create file objects from input parameters
reference_fp = file(params.reference)
Channel.fromFilePairs(params.reads, flat: true)
  .ifEmpty { exit 1, "error: did not find any read files value '${params.reads}'" }
  .set { ch_read_sets }


// Check reference exists
if (! reference_fp.exists()) {
  exit 1, "error: reference input '${reference_fp}' does not exist"
}


// TODO: check merge inputs exist and are complete i.e. no missing isolates
// TODO: ensure collisions in filename space between new and previous datasets
merge_source_dir = file(params.previous_run_dir)
merge_source_fastqc = Channel.fromPath(merge_source_dir / 'fastqc/individual_reports/*')
merge_source_gene_depth = Channel.fromPath(merge_source_dir / '*gene_depth.tsv')
merge_source_gene_coverage = Channel.fromPath(merge_source_dir / '*gene_coverage.tsv')
merge_source_mapping_stats = Channel.fromPath(merge_source_dir / '*mapping_stats.tsv')
merge_source_allele_matrices = Channel.fromPath(merge_source_dir / '*alleles.tsv')


// Run workflow - this unnamed workflow will be implicity executed
workflow {
  // TODO: investigate whether processes can access variables above scope
  // need to provide reference name for saving files to many processes - currently a little messy
  // params seem to work but only if they're set during start up i.e. in config/cmdline

  main:
    reference_data = prepare_reference(reference_fp)

    // Branch for SE/PE inputs here
    // NOTE: we can process both concurrently in one run - configuration would need to be very clear
    read_data = preprocess_reads(ch_read_sets, run_read_subsample, run_quality_assessment)
    variant_data = call_variants(read_data.reads, reference_data.fasta, reference_data.bt2_index, reference_data.samtools_index)

    // Mix SE/PE channels here
    // NOTE: just use the `mix` channel op
    stats_data = run_mapping_stats(variant_data, reference_fp)
    allele_matrix_data = run_allele_matrices(variant_data.bams, stats_data.snp_sites, stats_data.passing, reference_data.fasta)

    // Execute merge run if required, otherwise run post analysis as normal
    if (merge_run) {
      // TODO: this is a lot of argument verbiage, think about a better approach
      // We could `mix` here to half the number of arguments?
      merge_data = run_merge_results(read_data.fastqc, merge_source_fastqc,
                                     stats_data.gene_coverage, merge_source_gene_coverage,
                                     stats_data.gene_depth, merge_source_gene_depth,
                                     stats_data.stats, merge_source_mapping_stats,
                                     allele_matrix_data.matrices, merge_source_allele_matrices)
      //run_post(merge_data.fastqc, merge_data.allele_matrices, reference_fp, run_phylogeny)
    } else {
      run_post(read_data.fastqc, allele_matrix_data.matrices, reference_fp, run_phylogeny)
    }
}
