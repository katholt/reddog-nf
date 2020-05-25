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



include { prepare_reference;
          create_read_quality_reports;
          aggregate_read_quality_reports;
          subsample_reads;
          create_mpileups;
          calculate_gene_coverage_depth;
          calculate_mapping_statistics;
          aggregate_mapping_statistics;
          aggregate_snp_sites;
          create_allele_matrix;
          aggregate_allele_matrices;
          create_snp_alignment;
          determine_coding_consequences;
          infer_phylogeny } from './workflows/common.nf'

include { align_reads_to_reference;
          call_snps } from './workflows/paired_reads.nf'




workflow {
  // NOTE: these do not necessarily have to be nested in this scope
  // Create file objects from input parameters
  reference_fp = file("${params.reference}")

  // TODO: investigate whether processes can access variables above scope
  // need to provide reference name for saving files to many processes
  // params seem to work but only if they're set during start up i.e. in config/cmdline

  // TODO: work flow for paired-end run, single-end run, and merge run
  // conditionally branch here

  main:
    reference = prepare_reference(reference_fp)

    if (run_read_subsample) {
      ch_read_sets = subsample_reads(ch_read_sets)
    }

    if (run_quality_assessment) {
      fastqc_reports = create_read_quality_reports(ch_read_sets)
      aggregate_read_quality_reports(fastqc_reports.output.collect())
    }

    reads_aligned = align_reads_to_reference(ch_read_sets, reference.fasta, reference.bt2_index)
    mpileups = create_mpileups(reads_aligned.bams)

    bams_and_mpileups = reads_aligned.bams.join(mpileups.output)
    snps = call_snps(bams_and_mpileups, reference.fasta, reference.samtools_index)


    // TODO: place this into a function
    bams_vcfs_stats_and_metrics = reads_aligned.bams
      .join(snps.vcfs)
      .join(snps.coverage_depth)
      .join(reads_aligned.metrics)


    mapping_stats = calculate_mapping_statistics(bams_vcfs_stats_and_metrics)
    mapping_stats_aggregated = aggregate_mapping_statistics(mapping_stats.output.collect(), reference.name)

    snp_sites = aggregate_snp_sites(snps.sites.collect(), mapping_stats_aggregated.isolate_replicons)

    gene_coverage_depth = calculate_gene_coverage_depth(mpileups.output_noid.collect(), reference_fp)


    // TODO: place this into a function
    // Read in isolates the replicons that passed and generate a channel to emit [isolate_id, replicon_ids]
    // Where replicon_ids is a string with each replicon_id separated by a single space
    // Also get count of passing isolates so we can scale resource allocation if using SLURM executor
    isolate_replicons_passing = mapping_stats_aggregated.isolate_replicons.flatMap { filepath ->
        filepath.readLines().collect { line ->
          tokens = line.tokenize('\t')
          [tokens[0], tokens[1..-1].join(' ')]
        }
      }
    isolate_passing_count = isolate_replicons_passing.count()

    // Remove isolates that have no replicons that pass mapping criteria and then add SNP site file to each BAM
    // We perform this here so that we do not run jobs for isolates that have no passing replicons
    bams_and_sites = reads_aligned.bams.join(isolate_replicons_passing).combine(snp_sites.output)


    allele_matrices = create_allele_matrix(bams_and_sites, reference.fasta)


    // TODO: place this into a function
    // Create a channel that emits allele matrices arranged by replicon_id
    //   - input of [isolate_id, list(isolate_allele_matrices)]
    //     - nextflow returns a list for multiple files or single object for one file
    //     - check for different object types and process accordingly
    //   - use isolate_id to robustly get replicon_id from allele matrix filename
    //   - flat emit [replicon_id, isolate_allele_matrix] for each file
    //   - group each matrix by replicon_id to emit [replicon_id, list(isolate_allele_matrices)]
    replicon_allele_matrices = allele_matrices.output.flatMap { isolate_id, filepaths ->
        if (! (filepaths instanceof List)) {
          replicon_id = filepaths.getName().minus("_${isolate_id}_alleles.tsv")
          return [[replicon_id, filepaths]]
        } else {
          return filepaths.collect { filepath ->
              replicon_id = filepath.getName().minus("_${isolate_id}_alleles.tsv")
              [replicon_id, filepath]
            }
        }
      }.groupTuple()

    allele_matrices_aggregated = aggregate_allele_matrices(replicon_allele_matrices, snp_sites.output, reference.name)


    // TODO: place this into a function
    // Filter matrices that have no alleles so we don't needlessly execute downstream processes
    allele_matrices_filtered = allele_matrices_aggregated.output.filter { replicon_id, fp ->
        // Read first two lines of allele matrix and determine if we have data
        has_alleles = true
        fh = fp.newReader()
        for (int i = 0; i < 2; i++) { has_alleles = fh.readLine() != null }
        return has_alleles
      }

    determine_coding_consequences(allele_matrices_filtered, reference_fp)

    // TODO: restore isolate_count and run_phylogeny logic
    isolate_count = 500

    snp_alignment = create_snp_alignment(allele_matrices_filtered, reference.name)
    if (isolate_count <= 1000 | run_phylogeny) {
      infer_phylogeny(snp_alignment, reference.name)
    }
}
