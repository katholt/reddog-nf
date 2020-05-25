#!/usr/bin/env nextflow
// Enable DSL2
nextflow.preview.dsl=2


log.info """
‌‌

                              d8b           d8b
                              88P           88P
                             d88           d88
  88bd88b     d8888b     d888888       d888888       d8888b      d888b8b
  88P'  `    d8b_,dP    d8P' ?88      d8P' ?88      d8P' ?88    d8P' ?88
 d88         88b        88b  ,88b     88b  ,88b     88b  d88    88b  ,88b
d88'         `?888P'    `?88P'`88b    `?88P'`88b    `?8888P'    `?88P'`88b
                                                                       )88
                                                                      ,88P
                                                                  `?8888P

‌‌
""".stripIndent()


def print_help() {
  log.info """
  ==========================================================================
  reddog-nf: reddog nextflow implementation
  ==========================================================================
  Homepage (original pipeline): https://github.com/katholt/RedDog
  Homepage (nextflow implementation): https://github.com/scwatts/reddog-nf

  Required:
  --------------------
    --reads               Input reads (paired, gzip compressed)

    --reference           Reference genbank

    --output_dir          Output directory
  --------------------

  Other:
  --------------------
    --help                Displays this help message
  --------------------

  ==========================================================================

  """.stripIndent()
}


// Execute a command and capture standard streams along with return code
def execute_command(String command) {
  stdout = new StringBuilder()
  stderr = new StringBuilder()
  process = command.execute()
  process.waitForProcessOutput(stdout, stderr)
  return [process.exitValue(), stdout, stderr]
}


// For an optional stage param variable, check that it is either a Boolean or String
// If it is a string and either 'true' or 'false', return the boolean equivalent
def check_boolean_option(Object option, String name) {
  if (option.getClass() == java.lang.Boolean) {
    return option
  } else if (option.getClass() == java.lang.String) {
    if (option.toLowerCase() == 'true') {
      return true
    } else if (option.toLowerCase() == 'false') {
      return false
    }
  }
  exit 1, "error: ${name} option must be true or false"
}


// Check required input and outputs
if (params.help) {
  print_help()
  exit 0
}
if (! params.reads) {
  print_help()
  println "‌‌"
  exit 1, "error: option --reads is required"
}
if (! params.reference) {
  print_help()
  println "‌‌"
  exit 1, "error: option --reference is required"
}
if (! params.output_dir) {
  print_help()
  println "‌‌"
  exit 1, "error: option --output_dir is required"
}


// Require optional stage variables to be boolean
// We must check and change values if needed. The param variables are immutable so instead we declare new ones
run_read_subsample = check_boolean_option(params.subsample_reads, 'subsample_reads')
run_quality_assessment = check_boolean_option(params.quality_assessment, 'quality_assessment')
run_phylogeny = check_boolean_option(params.force_tree, 'force_tree')


// Check integer params
if (run_read_subsample & params.subsample_read_count.getClass() != java.lang.Integer) {
  if (! params.subsample_read_count.isInteger()) {
    exit 1, "error: subsample_read_count must be an integer, got {params.subsample_read_count}"
  }
}


// Create file objects from input parameters
reference_gbk_fp = file("${params.reference}")
// TODO: do we still need this?
output_dir = file("${params.output_dir}")


// Create channel for input read sets and get number of input isolates
Channel.fromFilePairs(params.reads, flat: true)
  .ifEmpty { exit 1, "error: did not find any read files value '${params.reads}'" }
  .set { ch_read_sets }


// Check reference exists
if (! reference_gbk_fp.exists()) {
  exit 1, "error: reference input '${reference_gbk_fp}' does not exist"
}


// Set path to project bin directory
bin_dir = "${workflow.projectDir}/bin"


// Validate reference
script = 'validate_reference.py'
validate_command = "${bin_dir}/${script} --reference_fp ${params.reference}"
(return_code, stdout, stderr) = execute_command(validate_command)
if (return_code != 0) {
  exit 1, "error: validation of reference failed and exited with the following message:\n\n${stderr}"
}


// TODO: add process and use local executor
/*
// Validate first 10 readsets
script = 'validate_reads.py'
ch_read_set_fps = ch_read_sets.take(10).flatMap { it -> it[1..2] }.collect()
validate_command = "${bin_dir}/${script} --reads_fps ${ch_read_set_fps.val.join(' ')}"
(return_code, stdout, stderr) = execute_command(validate_command)
if (return_code != 0) {
  exit 1, "error: validation of first ten reads failed and exited with the following message:\n\n${stderr}"
}
*/


// Do not run if output exists and contains files other than the run info directory (which is created by this point)
output_dir_files = []
output_dir.eachFile { output_dir_files.add(it.name) }
run_info_dirname = file(params.run_info_dir).simpleName
output_dir_files.remove(run_info_dirname)
if (output_dir_files.size() > 0 && ! params.force) {
  exit 1, "error: output directory '${output_dir}' already exists and contains other files, remove or use --force to overwrite"
}


// Do not run on MASSIVE unless user specifies profile to use to avoid inadvertently using a local executor
massive_hostnames = ['m3-login1', 'm3-login2']
on_massive = massive_hostnames.contains(InetAddress.getLocalHost().getHostName())
profile_explicit = workflow.commandLine.tokenize(' ').contains('-profile')
if (on_massive && ! profile_explicit) {
  exit 1, "error: to run on MASSIVE you must explicitly set -profile"
}


include preprocess_reads from './src/modules/read_preprocessing.nf'
include call_variants from './src/modules/variant_calling.nf'
include run_mapping_stats from './src/modules/mapping_stats.nf'
include run_allele_matrices from './src/modules/allele_matrices.nf'

include { prepare_reference; determine_coding_consequences; create_snp_alignment; infer_phylogeny } from './src/processes/common.nf'


workflow {
  // NOTE: these do not necessarily have to be nested in this scope
  // Create file objects from input parameters
  reference_fp = file("${params.reference}")

  // TODO: investigate whether processes can access variables above scope
  // need to provide reference name for saving files to many processes
  // params seem to work but only if they're set during start up i.e. in config/cmdline

  main:
    // TODO: work flow for paired-end run, single-end run, and merge run
    // conditionally branch here
    ch_read_sets = preprocess_reads(ch_read_sets, run_read_subsample, run_quality_assessment)
    reference_data = prepare_reference(reference_fp)
    variant_data = call_variants(ch_read_sets, reference_data.fasta, reference_data.bt2_index, reference_data.samtools_index)
    stats_data = run_mapping_stats(variant_data, reference_gbk_fp)
    allele_matrix_data = run_allele_matrices(variant_data.bams, stats_data.snp_sites,
                                              stats_data.isolate_replicons_passing, reference_data.fasta)

    // TODO: restore isolate_count and run_phylogeny logic
    isolate_count = 500

    snp_alignment = create_snp_alignment(allele_matrix_data.matrices, reference_data.name)
    if (isolate_count <= 1000 | run_phylogeny) {
      infer_phylogeny(snp_alignment, reference_data.name)
    }

    determine_coding_consequences(allele_matrix_data.matrices, reference_fp)
}
