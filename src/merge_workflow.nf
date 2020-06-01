// Processes
include gene_depth from './processes/merge.nf'
include gene_coverage from './processes/merge.nf'
include allele_matrix from './processes/merge.nf'
include mapping_stats from './processes/merge.nf'

// Channel helper functions
include get_replicon_id from './channel_helpers.nf'


workflow merge {
    take:
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
        if (run_quality_assessment) {
            // Symlink previous FastQC reports into new output directory
            fastqc_individual_output_dir = file(params.output_dir) / 'fastqc/individual_reports/'
            if (! fastqc_individual_output_dir.exists()) {
                fastqc_individual_output_dir.mkdirs()
            }
            merge_source_fastqc.map { filepath_src ->
                filepath_dst = fastqc_individual_output_dir / filepath_src.getName()
                if (! filepath_dst.exists()) {
                java.nio.file.Files.createSymbolicLink(filepath_dst, filepath_src)
                }
            }
            // Filter for zip files and add these to existing channel
            ch_fastqc = ch_fastqc.flatten().mix(merge_source_fastqc.filter { it.getName().endsWith('zip') })
        }

        // Merge gene stats tables
        gene_depth(gene_depth, merge_source_gene_depth, reference_name)
        gene_coverage(gene_coverage, merge_source_gene_coverage, reference_name)

        // Merge mapping stats
        merge_source_mapping_stats = get_replicon_id(merge_source_mapping_stats, '_mapping_stats.tsv', reference_name)
        mapping_stats = stats_data.map { [it.getName().minus('_mapping_stats.tsv'), it] }
        ch_mapping_stats_merge = mapping_stats.mix(merge_source_mapping_stats).groupTuple()
        mapping_stats(ch_mapping_stats_merge, reference_name)

        // TODO: catch where an allele matrix is create for a replicon in one run but not another
        //       just skip merging for these tables and send straight to filter channel
        // Merge allele matrices and update channel
        merge_source_allele_matrices = get_replicon_id(merge_source_allele_matrices, '_alleles.tsv', reference_name)
        ch_allele_matrices_merge = ch_allele_matrices.mix(merge_source_allele_matrices).groupTuple()
        allele_matrices = allele_matrix(ch_allele_matrices_merge, reference_name)

    emit:
        fastqc = ch_fastqc
        allele_matrices = allele_matrices
}
