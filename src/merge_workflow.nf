// Processes
include { copy_bams } from './processes/merge.nf'
include { copy_vcfs } from './processes/merge.nf'
include { copy_fastqc } from './processes/merge.nf'
include { index_bam } from './processes/merge.nf'
include { gene_depth } from './processes/merge.nf'
include { gene_coverage } from './processes/merge.nf'
include { mapping_stats } from './processes/merge.nf'
include { collect_snp_sites } from './processes/merge.nf'

// Channel helper functions
include { get_replicon_id } from './channel_helpers.nf'


workflow merge {
  take:
    merge_source_bams
    merge_source_vcfs
    ch_fastqc
    merge_source_fastqc
    gene_depth
    merge_source_gene_depth
    gene_coverage
    merge_source_gene_coverage
    stats_data
    merge_source_mapping_stats
    merge_source_allele_matrices
    run_read_quality_report
    reference_name

  main:
    // Copy BAMs and VCFs
    copy_vcfs(merge_source_vcfs)
    copy_bams(merge_source_bams)

    // Copy FastQC reports, then add zip output to existing channel
    if (run_read_quality_report) {
      copy_fastqc(merge_source_fastqc)
      ch_fastqc = ch_fastqc.flatten().mix(merge_source_fastqc.filter { it.getName().endsWith('zip') })
    }

    // Index BAMs
    bam_data = index_bam(merge_source_bams)

    // Merge gene stats tables
    gene_depth(gene_depth, merge_source_gene_depth, reference_name)
    gene_coverage(gene_coverage, merge_source_gene_coverage, reference_name)

    // Merge mapping stats
    merge_source_mapping_stats = get_replicon_id(merge_source_mapping_stats, '_mapping_stats.tsv', reference_name)
    mapping_stats = stats_data.flatten().map { [it.getName().minus('_mapping_stats.tsv'), it] }
    ch_mapping_stats_merge = mapping_stats.mix(merge_source_mapping_stats).groupTuple()
    // Sort files into channels based on whether they can be merged
    ch_mapping_stats_merge = ch_mapping_stats_merge.branch { replicon_id, fps ->
      merge: fps.size() == 2
      no_merge: fps.size() != 2
    }
    stats_data = mapping_stats(ch_mapping_stats_merge.merge, reference_name)
    stats_data.mix(ch_mapping_stats_merge.no_merge)

    // Collect SNP sites identified in previous run from the allele matrix
    merge_source_allele_matrices = get_replicon_id(merge_source_allele_matrices, '_alleles.csv', reference_name)
    merge_snp_sites = collect_snp_sites(merge_source_allele_matrices)

  emit:
    bams = bam_data
    fastqc = ch_fastqc
    snp_sites = merge_snp_sites
    mapping_stats = stats_data.output
}
