#!/usr/bin/env nextflow
// TODO: pre-flight checks
// TODO: decide approach for merge runs
// TODO: Provide config and allow options set from commandline


// File I/O
//params.reads = 'data/subset_reads/*_{1,2}_sub.fastq.gz'
params.reads = 'data/reads/*_{1,2}.fastq.gz'
params.reference = file('data/NCTC13753.gbk')
params.output_dir = file('output')


// Seed channels
Channel
  .fromFilePairs(params.reads, flat:true)
  .set { ch_read_sets }


// Prepare reference
// Convert to FASTA format
process reference_convert {
  output:
  file '*.fasta' into ch_reference_fasta

  script:
  reference_fasta = params.reference.simpleName + '.fasta'
  """
  genbank_to_fasta.py --input_fp ${params.reference} --output_fp ${reference_fasta}
  """
}
ch_reference_fasta.into {
    ch_index_reference_bt2;
    ch_index_reference_samtools;
    ch_collect_replicon_names;
    ch_align_reference;
    ch_snps_call_reference;
    ch_consensus_reference
  }


// Index reference for bowtie2
process reference_bowtie2_index {
  input:
  file reference from ch_index_reference_bt2

  output:
  file '*fasta.*bt2' into ch_reference_bt2_index

  script:
  """
  bowtie2-build ${reference} ${reference}
  """
}


// Index reference for samtools
process reference_samtools_index {
  input:
  file reference from ch_index_reference_samtools

  output:
  file '*fasta.fai' into ch_reference_samtools_index

  script:
  """
  samtools faidx ${reference}
  """
}


// Get replicons from FASTA reference we'll be mapping to
// TODO: python script get check replicon names then output to stdout
process collect_replicon_names {
  input:
  file reference from ch_collect_replicon_names

  output:
  env replicons into ch_replicons

  script:
  """
  replicons=\$(grep '^>' ${reference} | sed -e 's/^>//' -e 's/ .\\+\$//')
  """
}


// Align reads
process align_reads {
  input:
  tuple sample_id, file(reads_fwd), file(reads_rev) from ch_read_sets
  file reference from ch_align_reference
  file reference_indices from ch_reference_bt2_index

  output:
  file "${sample_id}.bam" into ch_bams

  script:
  """
  bowtie2 --sensitive-local -x ${reference} -1 ${reads_fwd} -2 ${reads_rev} -X 2000 | samtools view -ubS - | samtools sort - -o ${sample_id}.bam
  """
}
ch_bams.into { ch_sam_stats; ch_filter_unmapped }


// Get sam stats
// TODO: check where these are used
process calculate_sam_stats {
  input:
  file bam_fp from ch_sam_stats

  output:
  file '*_sam_stats.txt'

  script:
  sample_id = bam_fp.simpleName
  """
  sam-stats -A -B ${bam_fp} > ${sample_id}_sam_stats.txt
  """
}


// Filter unmapped reads from bam
process filter_unmapped_reads {
  publishDir "${params.output_dir}/bams/", saveAs: { fn -> fn.replaceFirst(/_filtered/, '') }

  input:
  file bam_fp from ch_filter_unmapped

  output:
  file '*_filtered.bam' into _ch_bams_filtered

  script:
  sample_id = bam_fp.simpleName
  """
  samtools view -hub -F 4 ${bam_fp} | samtools sort - -o ${sample_id}_filtered.bam
  """
}
_ch_bams_filtered.into { ch_index_filtered_bams; ch_bams_filtered }


// Index filtered bams
process index_filtered_bams {
  input:
  file bam_fp from ch_index_filtered_bams

  output:
  file '*bam.bai' into ch_bams_filtered_indices

  script:
  """
  samtools index ${bam_fp}
  """
}


// Create channel for bams and respective indices
ch_bams_filtered.merge(ch_bams_filtered_indices).into {
    ch_call_snps_bams_and_indices;
    ch_consensus_bams_and_indices;
    ch_coverage_bams_and_indices;
  }


// Call SNPs
// Create channel of combinations for replicion ~ isolate bam file
ch_replicons.tokenize(' ').flatMap().set { ch_call_snps_replicons }
ch_call_snps_replicons.combine(ch_call_snps_bams_and_indices).set { ch_call_snps }

process call_snps {
  input:
  tuple val(replicon), file(bam_fp), file(index_fp) from ch_call_snps
  file reference from ch_snps_call_reference
  file fasta_index from ch_reference_samtools_index

  output:
  file '*_raw.bcf' into ch_raw_bcfs

  script:
  sample_id = bam_fp.simpleName
  """
  samtools mpileup -u -t DP -f ${reference} ${bam_fp} -r ${replicon} | bcftools call -O b -cv - > ${sample_id}_raw.bcf
  """
}


// Get consensus sequences
process get_consensus {
  input:
  file reference from ch_consensus_reference
  tuple file(bam_fp), file(index_fp) from ch_consensus_bams_and_indices

  output:
  file '*cns.fq' into ch_consensus_sequences

  script:
  sample_id = bam_fp.simpleName
  """
  samtools mpileup -q 20 -ugB -f ${reference} ${bam_fp} | bcftools call -c - | vcfutils.pl vcf2fq > ${sample_id}_cns.fq
  """
}


// Get consensus sequences
process get_coverage {
  input:
  tuple file(bam_fp), file(index_fp) from ch_coverage_bams_and_indices

  output:
  file '*coverage.txt' into ch_coverage

  script:
  sample_id = bam_fp.simpleName
  """
  samtools mpileup ${bam_fp} | cut - -f 1-4 > ${sample_id}_coverage.txt
  """
}
