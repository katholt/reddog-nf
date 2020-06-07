// Processes
include gene_depth from './processes/merge.nf'
include gene_coverage from './processes/merge.nf'
include allele_matrix from './processes/merge.nf'
include mapping_stats from './processes/merge.nf'

// Channel helper functions
include get_replicon_id from './channel_helpers.nf'

// Utility functions
include symlink_merge_data from './utilities.nf'


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
        ch_allele_matrices
        merge_source_allele_matrices
        run_quality_assessment
        reference_name

    main:
        // Symlink BAMs and VCFs
        bam_output_dir = file(params.output_dir) / 'bams/'
        vcf_output_dir = file(params.output_dir) / 'vcfs/'
        symlink_merge_data(merge_source_bams, bam_output_dir)
        symlink_merge_data(merge_source_vcfs, vcf_output_dir)

        // Symlink FastQC reports, then add zip output to existing channel
        if (run_quality_assessment) {
            fastqc_individual_output_dir = file(params.output_dir) / 'fastqc/individual_reports/'
            symlink_merge_data(merge_source_fastqc, fastqc_individual_output_dir)
            ch_fastqc = ch_fastqc.flatten().mix(merge_source_fastqc.filter { it.getName().endsWith('zip') })
        }

        // Merge gene stats tables
        gene_depth(gene_depth, merge_source_gene_depth, reference_name)
        gene_coverage(gene_coverage, merge_source_gene_coverage, reference_name)

        // Merge mapping stats
        merge_source_mapping_stats = get_replicon_id(merge_source_mapping_stats, '_mapping_stats.tsv', reference_name)
        mapping_stats = stats_data.map { [it.getName().minus('_mapping_stats.tsv'), it] }
        ch_mapping_stats_merge = mapping_stats.mix(merge_source_mapping_stats).groupTuple()
        // Filter replicon stats which only had output in one of the two runs
        ch_mapping_stats_merge = ch_mapping_stats_merge.filter { replicon_id, fps -> fps.size() == 2 }
        mapping_stats(ch_mapping_stats_merge, reference_name)

        // Merge allele matrices and update channel
        merge_source_allele_matrices = get_replicon_id(merge_source_allele_matrices, '_alleles.tsv', reference_name)
        ch_allele_matrices_merge = ch_allele_matrices.mix(merge_source_allele_matrices).groupTuple()
        // Split replicon matrices by number of files found
        ch_allele_matrices_merge = ch_allele_matrices_merge.branch { replicon_id, fps ->
                multiple: fps.size() > 1
                single: fps.size() == 1
            }
        allele_matrices = allele_matrix(ch_allele_matrices_merge.multiple, reference_name)
        // Add non-merged matrices to output channel
        allele_matrices = allele_matrices.mix(ch_allele_matrices_merge.single)

    emit:
        fastqc = ch_fastqc
        allele_matrices = allele_matrices
}
