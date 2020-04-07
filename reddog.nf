#!/usr/bin/env nextflow
// TODO: with very large datasets we will error on too many command line arguments - one option will be to write to file in 'script' blocks

// TODO: pre-flight checks, see RedDog.py
// TODO: provide config and allow options set from commandline
// TODO: decide approach for merge runs

// TODO: sort inputs by size so the slowest jobs start first - this can provide small speed improvement

// TODO: create variant where we agglomerate stages in 3-4 steps which python scripts to handle execution


// File I/O
params.reads = 'data/reads/*_{1,2}.fastq.gz'
//params.reads = 'data/subset_reads/*_{1,2}_sub.fastq.gz'
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
    ch_consensus_reference;
    ch_genome_alignment_reference
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
process collect_replicon_names {
  input:
  file reference_fp from ch_collect_replicon_names

  output:
  env replicons into ch_replicons

  script:
  """
  # Get replicon names then check they are all unique
  replicons=\$(grep '^>' ${reference_fp} | sed -e 's/^>//' -e 's/ .\\+\$//')
  replicons_duplicate=\$(echo "\${replicons}" | sort | uniq -d)
  replicons_duplicate_count=\$(echo -n "\${replicons_duplicate}" | sort | uniq -d | wc -l)
  if [[ "\${replicons_duplicate_count}" -gt 0 ]]; then
    echo "Found duplicate replicon names, you must rename these" 1>&2
    echo "Duplicates: \${replicons_duplicate}" 1>&2
    exit 1
  fi
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
  tuple val(sample_id), file("${sample_id}.bam") into ch_raw_bams

  script:
  """
  bowtie2 --sensitive-local -x ${reference_fp} -1 ${reads_fwd} -2 ${reads_rev} -X 2000 | samtools view -ubS - | samtools sort - -o ${sample_id}.bam
  """
}
ch_raw_bams.into { ch_filter_unmapped; ch_replicon_stats_bams }


// Filter unmapped reads from bam
process filter_unmapped_reads {
  publishDir "${params.output_dir}/bams/"

  input:
  tuple sample_id, file(bam_fp) from ch_filter_unmapped

  output:
  tuple val(sample_id), file('*_filtered.bam') into ch_bams_filtered

  script:
  """
  samtools view -hub -F 4 ${bam_fp} | samtools sort - -o ${sample_id}_filtered.bam
  """
}


// Index filtered bams
process index_filtered_bams {
  input:
  tuple sample_id, file(bam_fp) from ch_bams_filtered

  output:
  tuple val(sample_id), file(bam_fp), file('*bam.bai') into ch_bams_and_indices

  script:
  """
  samtools index ${bam_fp}
  """
}
ch_bams_and_indices.into {
    ch_call_snps_bams_and_indices;
    ch_consensus_bams_and_indices
  }


// Call SNPs
// Create channel of combinations for replicion ~ isolate bam files
ch_call_snps_replicons.combine(ch_call_snps_bams_and_indices).set { ch_call_snps }
process call_snps {
  input:
  tuple replicon_id, sample_id, file(bam_fp), file(index_fp) from ch_call_snps
  file reference_fp from ch_snps_call_reference
  file fasta_index from ch_reference_samtools_index

  output:
  tuple val(sample_id), val(replicon_id), file('*_raw.bcf') into ch_raw_bcfs

  script:
  """
  samtools mpileup -u -t DP -f ${reference_fp} ${bam_fp} -r ${replicon_id} | bcftools call -O b -cv - > ${sample_id}_${replicon_id}_raw.bcf
  """
}


// Call consensus
process call_consensus {
  input:
  tuple sample_id, file(bam_fp), file(index_fp) from ch_consensus_bams_and_indices
  file reference_fp from ch_consensus_reference

  output:
  file '*consensus.fastq' into ch_consensus

  script:
  """
  samtools mpileup -q 20 -ugB -f ${reference_fp} ${bam_fp} | bcftools call -c - | vcfutils.pl vcf2fq > "${sample_id}_consensus.fastq"
  """
}
// We need a single list item emitted to combine with VCFs for the allele matrix creatation step
// This works but there is likely a better way
ch_consensus.toList().toList().set { ch_allele_matrix_consensus }


// Filter SNPs using vcfutils
process filter_snps_vcfutils {
  input:
  tuple sample_id, replicon_id, file(bcf_fp) from ch_raw_bcfs

  output:
  tuple val(sample_id), val(replicon_id), file('*_filtered.vcf') into ch_filter_snps_q30_hets

  script:
  """
  bcftools view ${bcf_fp} | vcfutils.pl varFilter -d 5 -D 23 -Q 30 > ${sample_id}_${replicon_id}_filtered.vcf
  """
}


// Filter SNPs with less than Q30 and separate heterozygous SNPs
process filter_snps_q30_hets {
  publishDir "${params.output_dir}/vcfs/"

  input:
  tuple sample_id, replicon_id, file(vcf_fp) from ch_filter_snps_q30_hets

  output:
  tuple val(sample_id), val(replicon_id), file('*_q30.vcf'), file('*hets.vcf') into ch_filtered_vcfs

  script:
  """
  filter_snp_calls.py --input_vcf_fp ${vcf_fp} --output_q30_vcf_fp ${sample_id}_${replicon_id}_q30.vcf --output_hets_vcf_fp ${sample_id}_${replicon_id}_hets.vcf
  """
}
// Create channels that emit sample_id, replicon_id, q30_vcf_fp, hets_vcf_fp (ordered by replicon_id)
// NOTE: replicon ordering is required here so that the calculate_replicon_statistics process output glob
// matches input. This also means that sort order resulting from toSortedList MUST match glob sort order.
// An assertion is provided in the aggregation process
// These channels emit sample_id, list(replicon_ids), list(q30), list(hets)
ch_filtered_vcfs.toSortedList { items -> items[1] }.flatMap().into { ch_replicon_stats_vcfs; _ch_allele_matrix_vcfs }


// Calculate replicon statistics
// Group VCFs by sample_id and add respective BAM filepath
ch_replicon_stats_vcfs.groupTuple().join(ch_replicon_stats_bams).set { ch_replicon_stats }
process calculate_replicon_statistics {
  input:
  tuple sample_id, replicons, file(vcf_q30_fps), file(vcf_hets_fps), file(bam_fp) from ch_replicon_stats

  output:
  tuple val(replicons), file("*${replicons}_RepStats.txt") into ch_replicon_stats_per_sample

  script:
  """
  calculate_replicon_stats.py --raw_bam_fp ${bam_fp} --vcf_q30_fps ${vcf_q30_fps} --vcf_hets_fps ${vcf_hets_fps} --output_dir ./
  """
}


// Aggregate replicon statistics
// Transpost channel such that we group stats files of the same replicon
ch_replicon_stats_per_sample.transpose().groupTuple().set { ch_replicon_stats_aggregate }
process aggregate_replicon_statistics {
  publishDir "${params.output_dir}", saveAs: { filename -> "${params.reference.simpleName}_${filename}" }

  input:
  tuple replicon_id, file(replicon_stats_fps) from ch_replicon_stats_aggregate

  output:
  file '*RepStats.txt'

  script:
  """
  # First check we have the corret replicon stats file
  replicon_id_found=\$(head -n1 -q ${replicon_stats_fps} | sed 's/^#//' | sort | uniq)
  replicon_id_counts=\$(echo "\${replicon_id_found}" | wc -l)
  if [[ "\${replicon_id_counts}" -gt 1 ]]; then
    echo 'error: got more than one replicon during aggregation' 1>&2
    exit 1;
  elif [[ "\${replicon_id_found}" != "${replicon_id}" ]]; then
    echo 'error: the wrong replicon during aggregation - expected ${replicon_id}, got \${replicon_id_found}' 1>&2
    exit 1;
  fi;

  sed -n '2p' ${replicon_stats_fps} > ${replicon_id}_RepStats.txt
  sed -e '/^#/d' -e '/^Isolate/d' ${replicon_stats_fps} | get_ingroup_outgroup.awk >> ${replicon_id}_RepStats.txt
  """
}


// Create allele matrix for each replicon
// Transform channel to only emit replicon_id and q30 VCFs (grouped by replicon)
_ch_allele_matrix_vcfs.groupTuple(by: 1).map { items -> items[1..2] }.set { ch_allele_matrix_vcfs }
ch_allele_matrix_vcfs.combine(ch_allele_matrix_consensus).set { ch_create_allele_matrix }
process create_allele_matrix {
  publishDir "${params.output_dir}", saveAs: { filename -> "${params.reference.simpleName}_${filename}" }

  input:
  tuple replicon_id, file(vcf_fps), file(consensus_fps) from ch_create_allele_matrix

  output:
  tuple val(replicon_id), file('*_alleles_var.tsv') into ch_filter_allele_matrix

  script:
  """
  create_allele_matrix.py --vcf_fps ${vcf_fps} --consensus_fps ${consensus_fps} --replicon ${replicon_id} > ${replicon_id}_alleles_var.tsv
  """
}

// Filter allele matrix
process filter_allele_matrix {
  publishDir "${params.output_dir}", saveAs: { filename -> "${params.reference.simpleName}_${filename}" }

  input:
  tuple replicon_id, file(allele_matrix_fp) from ch_filter_allele_matrix

  output:
  tuple val(replicon_id), file("*_cons0.95.tsv") into ch_filtered_allele_matrix

  script:
  """
  filter_allele_matrix.py --allele_fp ${allele_matrix_fp} > ${replicon_id}_alleles_var_cons0.95.tsv
  """
}
ch_filtered_allele_matrix.into { ch_genome_alignment; ch_determine_coding_consequences }


// Create pseudo genome alignment
process create_pseudo_genome_alignment {
  publishDir "${params.output_dir}", saveAs: { filename -> "${params.reference.simpleName}_${filename}" }

  input:
  tuple replicon_id, file(allele_matrix_fp) from ch_genome_alignment
  file reference_fp from ch_genome_alignment_reference

  output:
  tuple val(replicon_id), file("*.mfasta") into ch_infer_phylogeny

  script:
  """
  create_pseudo_genome_alignment.py --allele_fp ${allele_matrix_fp} --reference_fp ${reference_fp} --replicon ${replicon_id} > ${replicon_id}_alleles_var_cons0.95.mfasta
  """
}


// Determine coding consequences
process determine_coding_consequences {
  publishDir "${params.output_dir}", saveAs: { filename -> "${params.reference.simpleName}_${filename}" }

  input:
  tuple replicon_id, file(allele_matrix_fp) from ch_determine_coding_consequences

  output:
  file "*consequences.tsv"

  script:
  """
  determine_coding_consequences.py --allele_fp ${allele_matrix_fp} --reference ${params.reference} --replicon ${replicon_id} > ${replicon_id}_alleles_var_cons0.95_consequences.tsv
  """
}


// Infer phylogeny
process infer_phylogeny {
  publishDir "${params.output_dir}", saveAs: { filename -> "${params.reference.simpleName}_${filename}" }

  input:
  tuple replicon_id, file(alignment_fp) from ch_infer_phylogeny

  output:
  file '*.tree'

  script:
  """
  FastTree -gtr -gamma -nt ${alignment_fp} > ${replicon_id}_alleles_var_cons0.95.tree
  """
}
