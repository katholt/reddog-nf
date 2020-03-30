#!/usr/bin/env nextflow
// TODO: pre-flight checks
// TODO: decide approach for merge runs
// TODO: provide config and allow options set from commandline


// File I/O
//params.reads = 'data/reads/*_{1,2}.fastq.gz'
params.reads = 'data/subset_reads/*_{1,2}_sub.fastq.gz'
params.reference = file('data/NCTC13753.gbk')
params.output_dir = file('output')


// Create channel for input read sets
Channel.fromFilePairs(params.reads, flat:true).set { ch_read_sets }


// Convert reference to FASTA format
process convert_reference {
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
process index_reference_bowtie2 {
  input:
  file reference_fp from ch_index_reference_bt2

  output:
  file '*fasta.*bt2' into ch_reference_bt2_index

  script:
  """
  bowtie2-build ${reference_fp} ${reference_fp}
  """
}


// Index reference for samtools
process index_reference_samtools {
  input:
  file reference_fp from ch_index_reference_samtools

  output:
  file '*fasta.fai' into ch_reference_samtools_index

  script:
  """
  samtools faidx ${reference_fp}
  """
}


// Get replicons from FASTA reference we'll be mapping to
// TODO: python script get check replicon names then output to stdout
process collect_replicon_names {
  input:
  file reference_fp from ch_collect_replicon_names

  output:
  env replicons into ch_replicons

  script:
  """
  replicons=\$(grep '^>' ${reference_fp} | sed -e 's/^>//' -e 's/ .\\+\$//')
  """
}
ch_replicons.tokenize(' ').flatMap().set { ch_call_snps_replicons }


// Align reads
process align_reads_to_reference {
  input:
  tuple sample_id, file(reads_fwd), file(reads_rev) from ch_read_sets
  file reference_fp from ch_align_reference
  file reference_indices from ch_reference_bt2_index

  output:
  file "${sample_id}.bam" into ch_bams

  script:
  """
  bowtie2 --sensitive-local -x ${reference_fp} -1 ${reads_fwd} -2 ${reads_rev} -X 2000 | samtools view -ubS - | samtools sort - -o ${sample_id}.bam
  """
}
ch_bams.into { ch_sam_stats; ch_filter_unmapped }


// Filter unmapped reads from bam
process filter_unmapped_reads {
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
_ch_bams_filtered.into { ch_bams_filtered; ch_index_filtered_bams }


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


// Create channel to emit bams with the respective index
ch_bams_filtered.merge(ch_bams_filtered_indices).into {
    ch_call_snps_bams_and_indices;
    ch_consensus_bams_and_indices;
    ch_depth_coverage_bams_and_indices;
  }


// Call SNPs
// Create channel of combinations for replicion ~ isolate bam files
ch_call_snps_replicons.combine(ch_call_snps_bams_and_indices).set { ch_call_snps }
process call_snps {
  input:
  tuple val(replicon), file(bam_fp), file(index_fp) from ch_call_snps
  file reference_fp from ch_snps_call_reference
  file fasta_index from ch_reference_samtools_index

  output:
  file '*_raw.bcf' into ch_raw_bcfs

  script:
  sample_id = bam_fp.simpleName
  """
  samtools mpileup -u -t DP -f ${reference_fp} ${bam_fp} -r ${replicon} | bcftools call -O b -cv - > ${sample_id}_raw.bcf
  """
}


// Filter SNPs using vcfutils
process filter_snps_vcfutils {
  input:
  file bcf_fp from ch_raw_bcfs

  output:
  file '*_q30.vcf' into ch_filter_snps_q30_hets

  script:
  sample_id = bcf_fp.simpleName
  """
  bcftools view ${bcf_fp} | vcfutils.pl varFilter -d 5 -D 23 -Q 30 > ${sample_id}_q30.vcf
  """
}


// Filter SNPs with less than Q30 and separate heterozygous SNPs
process filter_snps_q30_hets {
  input:
  file vcf_fp from ch_filter_snps_q30_hets

  output:
  file '*_q30.vcf' into ch_filtered_vcfs

  script:
  """
  filter_snp_calls.py --input_vcf_fp ${vcf_fp} --output_q30_vcf_fp ${sample_id}_q30.vcf --output_hets_vcf_fp ${sample_id}_hets.vcf
  """
}


// Get consensus sequences
process get_consensus {
  input:
  file reference_fp from ch_consensus_reference
  tuple file(bam_fp), file(index_fp) from ch_consensus_bams_and_indices

  output:
  file '*cns.fq' into ch_consensus_sequences

  script:
  sample_id = bam_fp.simpleName
  """
  samtools mpileup -q 20 -ugB -f ${reference_fp} ${bam_fp} | bcftools call -c - | vcfutils.pl vcf2fq > ${sample_id}_cns.fq
  """
}


/*
// Get mapping depth and coverage
// TODO: remove this - we'll to it live (in python)
process get_mapping_depth_coverage {
  input:
  tuple file(bam_fp), file(index_fp) from ch_depth_coverage_bams_and_indices

  output:
  file '*depth_coverage.tsv' into ch_depth_coverage

  script:
  sample_id = bam_fp.simpleName
  """
  samtools mpileup -aa NCTC13753_*.bam | cut -f1,4 -d\$'\t' | \
    awk 'BEGIN {
           OFS="\t"
         }
         {
           sums[\$1] += \$2
           size[\$1] += 1
           if (\$2 > 0) {
             bases[\$1] += 1
           }
         } END {
           print "replicon", "replicon_size", "average_depth", "coverage"
           for (replicon in sums) {
             if (bases[replicon] == 0) {
               average_depth = 0.0
             } else {
               average_depth = sums[replicon] / bases[replicon]
             }
             print replicon, size[replicon], average_depth, bases[replicon] / size[replicon] * 100
           }
         }' > ${sample_id}_depth_coverage.tsv
  """
}
*/


/*
process get_vcf_stats {
  script:
  """
  awk 'BEGIN {
    OFS="\t"
    snps = twoalts = indels = 0
  } $1 !~ /^#/ {
    if ($8 ~ /^I/) {
      indels += 1
    } else {
      snps += 1
      if ($4 ~ /,/) {
        twoalts += 1
      }
    }
  } END {
    print "snps", "twoalts", "indels"
    print snps, twoalts, indels
  }' ${vcf_fp} > ${sample_id}_vcf_stats.tsv
  """
}
*/


/*
// Get sam stats
// TODO: remove this and associated dependency - we'll to it live (in python)
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
*/


/*
// Infer phylogeny
process infer_phylogeny {
  input:
  file alignment_fp from ch_alignments

  output:
  file '*_phylogeny.tree'

  script:
  replicon = alignment_fp.simpleName
  """
  FastTree -gtr -gamma -nt ${alignment_fp} > ${replicon}_phylogeny.tree
  """
}
*/
