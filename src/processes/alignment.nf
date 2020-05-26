// Align PE reads and apply some post-processing
//   - align with bowtie2
//     - minimum fragment length set to 2000, default of 500 too restrictive
//     - write metrics to output file, limit reporting to 30 minute intervals and final summary
//   - filtered unmapped reads
//   - sort reads
//   - index resulting BAM file
process align_reads_pe {
  publishDir "${params.output_dir}/bams/", pattern: '*.bam'

  input:
  tuple isolate_id, path(reads_fwd), path(reads_rev)
  path reference_fp
  path reference_indices

  output:
  tuple val(isolate_id), path("${isolate_id}.bam"), path("${isolate_id}.bam.*"), emit: bams
  tuple val(isolate_id), path("${isolate_id}_metrics.tsv"), emit: metrics

  script:
  """
  bowtie2 --sensitive-local -x ${reference_fp} -1 ${reads_fwd} -2 ${reads_rev} -X 2000 --met 1800 --met-file ${isolate_id}_metrics.tsv | \
    samtools view -u -F4 | samtools sort > ${isolate_id}.bam
  samtools index ${isolate_id}.bam
  """
}


// Align SE reads and apply some post-processing
//   - align with bowtie2
//     - minimum fragment length set to 2000, default of 500 too restrictive
//     - write metrics to output file, limit reporting to 30 minute intervals and final summary
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
  tuple val(isolate_id), path("${isolate_id}_metrics.tsv"), emit: metrics

  script:
  """
  bowtie2 --sensitive-local -x ${reference_fp} -U ${reads_fp} -X 2000 --met 1800 --met-file ${isolate_id}_metrics.tsv | \
    samtools view -u -F4 | samtools sort > ${isolate_id}.bam
  samtools index ${isolate_id}.bam
  """
}
