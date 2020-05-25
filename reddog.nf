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



// TODO: restore preflight checks



include { prepare_reference;
          create_mpileups;
          calculate_gene_coverage_depth;
          calculate_mapping_statistics;
          aggregate_mapping_statistics;
          aggregate_snp_sites;
          create_allele_matrix;
          aggregate_allele_matrices;
          create_snp_alignment;
          determine_coding_consequences;
          infer_phylogeny } from './src/common.nf'

include { align_reads_to_reference;
          call_snps } from './src/paired_reads.nf'




workflow {
  // NOTE: these do not necessarily have to be nested in this scope
  // Create file objects from input parameters
  reference_fp = file("${params.reference}")

  // TODO: investigate whether processes can access variables above scope
  // need to provide reference name for saving files to many processes
  // params seem to work but only if they're set during start up i.e. in config/cmdline

  // TODO: work flow for paired-end run, single-end run, and merge run
  // conditionally branch here

  // Create channel for input read sets and get number of input isolates
  Channel.fromFilePairs(params.reads, flat: true).set { ch_read_sets }

  main:
    reference = prepare_reference(reference_fp)
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
    run_phylogeny = true

    snp_alignment = create_snp_alignment(allele_matrices_filtered, reference.name)
    if (isolate_count <= 1000 | run_phylogeny) {
      infer_phylogeny(snp_alignment, reference.name)
    }
}
