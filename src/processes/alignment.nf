// Align PE reads and apply some post-processing
//   - align with bowtie2
//     - minimum fragment length set to 2000, default of 500 too restrictive
//   - split reads into mapped or unmapped BAMs
//     - sort reads in the mapped reads BAM during split
//   - index mapped reads BAM
process align_reads_pe {
  publishDir "${params.output_dir}/bams/", pattern: '*.bam'

  input:
  tuple isolate_id, path(reads_fwd), path(reads_rev)
  path reference_fp
  path reference_indices

  output:
  tuple val(isolate_id), path("${isolate_id}.bam"), path("${isolate_id}.bam.*"), emit: bams
  tuple val(isolate_id), path("${isolate_id}_unmapped.bam"), emit: bams_unmapped

  script:
  """
  bowtie2 ${params.bt2_mode} -x ${reference_fp} -1 ${reads_fwd} -2 ${reads_rev} -X ${params.bt2_max_frag_len} | \
    split_mapped_unmapped_reads.awk -v mapped_fp=${isolate_id}.bam -v unmapped_fp=${isolate_id}_unmapped.bam
  samtools index ${isolate_id}.bam
  """
}


// Align SE reads and apply some post-processing
//   - align with bowtie2
//     - minimum fragment length set to 2000, default of 500 too restrictive
//   - filtered unmapped reads
//   - sort reads
//   - index resulting BAM file
process align_reads_se {
  publishDir "${params.output_dir}/bams/", pattern: '*.bam'

  input:
  tuple isolate_id, path(reads_fp)
  path reference_fp
  path reference_indices

  output:
  tuple val(isolate_id), path("${isolate_id}.bam"), path("${isolate_id}.bam.*"), emit: bams
  tuple val(isolate_id), path("${isolate_id}_unmapped.bam"), emit: bams_unmapped

  script:
  """
  bowtie2 ${params.bt2_mode} -x ${reference_fp} -U ${reads_fp} -X ${params.bt2_max_frag_len} | \
    split_mapped_unmapped_reads.awk -v mapped_fp=${isolate_id}.bam -v unmapped_fp=${isolate_id}_unmapped.bam
  samtools index ${isolate_id}.bam
  """
}
