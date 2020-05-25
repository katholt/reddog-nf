#!/usr/bin/env nextflow
// Enable DSL2
nextflow.preview.dsl=2


// Workflows
include preprocess_reads from './src/modules/read_preprocessing.nf'
include call_variants from './src/modules/variant_calling.nf'
include run_mapping_stats from './src/modules/mapping_stats.nf'
include run_allele_matrices from './src/modules/allele_matrices.nf'

// Processes used in main workflow - separated for readability
include prepare_reference from './src/processes/common.nf'
include determine_coding_consequences from './src/processes/common.nf'
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
// We must check and change values if needed. The param variables are immutable so instead we declare new ones
run_read_subsample = check_boolean_option(params.subsample_reads, 'subsample_reads')
run_quality_assessment = check_boolean_option(params.quality_assessment, 'quality_assessment')
run_phylogeny = check_boolean_option(params.force_tree, 'force_tree')
merge_run = check_boolean_option(params.force_tree, 'merge_run')


// Validate merge run options and files in previous run
// TODO: implement this


// Check integer params
if (run_read_subsample & params.subsample_read_count.getClass() != java.lang.Integer) {
  if (! params.subsample_read_count.isInteger()) {
    exit 1, "error: subsample_read_count must be an integer, got {params.subsample_read_count}"
  }
}


// Create file objects from input parameters
reference_fp = file("${params.reference}")
Channel.fromFilePairs(params.reads, flat: true)
  .ifEmpty { exit 1, "error: did not find any read files value '${params.reads}'" }
  .set { ch_read_sets }


// Check reference exists
if (! reference_fp.exists()) {
  exit 1, "error: reference input '${reference_fp}' does not exist"
}


// Run workflow - this unnamed will be implicity executed
workflow {
  // NOTE: this do not necessarily have to be nested in this scope
  // Create file objects from input parameters

  // TODO: investigate whether processes can access variables above scope
  // need to provide reference name for saving files to many processes
  // params seem to work but only if they're set during start up i.e. in config/cmdline

  main:
    reference_data = prepare_reference(reference_fp)

    // Branch for SE/PE inputs
    // We can process both in a single run - configuration would need to be very clear
    ch_read_sets = preprocess_reads(ch_read_sets, run_read_subsample, run_quality_assessment)
    variant_data = call_variants(ch_read_sets, reference_data.fasta, reference_data.bt2_index, reference_data.samtools_index)

    // Mix SE/PE channels

    // TODO: update aggregate_allele_matrices to decouple core allele filtering
    stats_data = run_mapping_stats(variant_data, reference_fp)
    allele_matrix_data = run_allele_matrices(variant_data.bams, stats_data.snp_sites, stats_data.passing, reference_data.fasta)

    merge_run = false
    if (merge_run) {
      // Tasks:
      //   - read quality reports
      //     - determine if both runs have reports
      //     - run multiqc again
      //   - gene coverage, depth
      //     - merge only
      //   - mapping stats
      //     - merge and recompute ingroup/outgroup
      //   - allele matrices
      //     - merge only (filter later)

      // For now all data from previous run will be copied and source directory untouched
      //   - user to manually delete
    }

    // TODO: process here to filter core SNPs

    determine_coding_consequences(allele_matrix_data.matrices, reference_fp)

    // TODO: restore isolate_count logic - issues with channel->int in DSL2
    isolate_count = 500
    snp_alignment = create_snp_alignment(allele_matrix_data.matrices, reference_data.name)
    if (isolate_count <= 1000 | run_phylogeny) {
      infer_phylogeny(snp_alignment, reference_data.name)
    }
}
