#!/usr/bin/env nextflow
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


// Check arguments
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


// Create file objects from input parameters
reference = file("${params.reference}")
output_dir = file("${params.output_dir}")


// Create channel for input read sets
Channel.fromFilePairs(params.reads, flat:true).set { ch_read_sets }


// Convert reference to FASTA format
process convert_reference {
  output:
  file '*.fasta' into ch_reference_fasta

  script:
  reference_fasta = reference.simpleName + '.fasta'
  """
  genbank_to_fasta.py --input_fp ${reference} --output_fp ${reference_fasta}
  """
}
ch_reference_fasta.into {
    ch_index_reference_bt2;
    ch_index_reference_samtools;
    ch_collect_replicon_names;
    ch_align_reference;
    ch_call_snps_reference;
    ch_allele_matrix_reference;
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
  publishDir "${output_dir}/bams/"

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
    ch_allele_matrix_bams_and_indices
  }


// Call SNPs
// Create channel of combinations for replicion ~ isolate bam files
ch_call_snps_replicons.combine(ch_call_snps_bams_and_indices).set { ch_call_snps }

process call_snps {
  input:
  tuple replicon_id, sample_id, file(bam_fp), file(index_fp) from ch_call_snps
  file reference_fp from ch_call_snps_reference
  file fasta_index from ch_reference_samtools_index

  output:
  tuple val(sample_id), val(replicon_id), file('*_raw.bcf') into ch_raw_bcfs

  script:
  """
  samtools mpileup -u -t DP -f ${reference_fp} ${bam_fp} -r ${replicon_id} | bcftools call -O b -cv - > ${sample_id}_${replicon_id}_raw.bcf
  """
}


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
  publishDir "${output_dir}/vcfs/"

  input:
  tuple sample_id, replicon_id, file(vcf_fp) from ch_filter_snps_q30_hets

  output:
  tuple val(sample_id), val(replicon_id), file('*_q30.vcf'), file('*hets.vcf') into _ch_replicon_stats_vcfs
  tuple val(sample_id), val(replicon_id), file('*_q30.vcf') into ch_snp_sites_vcfs

  script:
  """
  filter_snp_calls.py --input_vcf_fp ${vcf_fp} --output_q30_vcf_fp ${sample_id}_${replicon_id}_q30.vcf --output_hets_vcf_fp ${sample_id}_${replicon_id}_hets.vcf
  """
}
// This channel emits sample_id, list(replicon_ids), list(q30), list(hets)
// NOTE: replicon ordering is required here so that the calculate_replicon_statistics process output glob
// matches input. This also means that sort order resulting from toSortedList MUST match glob sort order.
// An assertion is provided in the aggregation process
_ch_replicon_stats_vcfs.toSortedList { items -> items[1] }.flatMap().set{ ch_replicon_stats_vcfs }


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
  publishDir "${output_dir}", saveAs: { filename -> "${reference.simpleName}_${filename}" }

  input:
  tuple replicon_id, file(replicon_stats_fps) from ch_replicon_stats_aggregate

  output:
  file('*RepStats.txt') into ch_allele_matrix_rep_stats
  env sample_status into _ch_samples_status

  script:
  """
  # First check we have the correct replicon stats file
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
  sample_status=\$(awk -v var="${replicon_id}" 'NR != 1 { print var, \$1, \$10 }' ${replicon_id}_RepStats.txt)
  """
}
// This creates a channel that emits sample_id, replicon_id, pass/fail
// TODO: A dicey approach as it relies on every entry producing 3 data points, make more robust
_ch_samples_status.map { items -> items.tokenize(' ') }.flatMap().buffer(size: 3).set { ch_samples_status }


// Get SNP sites for each sample
process collect_snp_sites {
  input:
  tuple sample_id, replicon_id, file(vcf_fp) from ch_snp_sites_vcfs

  output:
  tuple val(replicon_id), file('*_sites.tsv') into ch_snp_sites_positions

  script:
  """
  get_snp_sites.awk ${vcf_fp} > ${replicon_id}_${sample_id}_sites.tsv
  """
}
// Group SNP sites by replicon id
ch_snp_sites_positions.groupTuple().set { ch_snp_sites_aggregate }


// Combine SNP sites
process aggregate_snp_sites {
  input:
  tuple replicon_id, file(snp_sites_fps) from ch_snp_sites_aggregate

  output:
  tuple val(replicon_id), file('*_sites.tsv') into ch_snp_sites

  script:
  """
  sort --merge -n --parallel=1 ${snp_sites_fps} | uniq > ${replicon_id}_sites.tsv
  """
}


// Here we aim to get a list of input files while excluding construction of allele matrix for failed samples
// This must be done on a per replicon basis - i.e. a sample may fail on one replicon but not others
// To achieve this only samples that pass mapping for a replicon are assigned the respective replicon SNP site
// file. This works because the script that constructs the allele matrix creates them from SNP site files; if
// a SNP site file is missing for a replicon then that matrix is not constructed.
// Channel emits sample_id, list(snp_site_fps) - for passing samples and replicons only. Process:
//   - remove all samples that failed to pass mapping criteria for each replicon
//   - samples passing for each sample are assigned a SNP site file
//   - remove replicon_id and then group by sample_id
ch_samples_status.filter { items -> items[2] != 'f' }
  .map { items -> items[0..1] }
  .combine(ch_snp_sites, by: 0)
  .map { items -> items[1,2] }
  .groupTuple()
  .set { ch_snp_sites_pass }

// Now add respective BAM and index to each set of SNP site files
ch_allele_matrix_bams_and_indices.join(ch_snp_sites_pass, by: 0).set { ch_create_allele_matrix }

process create_allele_matrix {
  input:
  tuple sample_id, file(bam_fp), file(index_fp), file(sites_fp) from ch_create_allele_matrix
  file reference_fp from ch_allele_matrix_reference

  output:
  tuple val(sample_id), file('*_alleles.tsv') into _ch_allele_matrix_aggregate

  script:
  """
  create_allele_matrix.py --bam_fp ${bam_fp} --sites_fps ${sites_fp} --reference_fp ${reference_fp} --output_dir ./
  """
}
_ch_allele_matrix_aggregate.println()


// TODO: get output of create_allele_matrix to match replicon_id to group for aggregation
// Get be robust (and slightly less than optimal) and grab the replicon_id from the file name

// TODO: Output is going to be a string for single items or list for more ugh
// We need to learn more groovy to use in closures for this


/*
// Aggregate allele matrices
process aggregate_allele_matrices {
  publishDir "${output_dir}", saveAs: { filename -> "${reference.simpleName}_${filename}" }

  input:
  file(allele_fps) from ch_allele_matrix_aggregate

  output:
  file('*alleles.tsv')

  script:
  """
  # some BASH
  """
}

/*

// Filter allele matrix
process filter_allele_matrix {
  publishDir "${output_dir}", saveAs: { filename -> "${reference.simpleName}_${filename}" }

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
  publishDir "${output_dir}", saveAs: { filename -> "${reference.simpleName}_${filename}" }

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
  publishDir "${output_dir}", saveAs: { filename -> "${reference.simpleName}_${filename}" }

  input:
  tuple replicon_id, file(allele_matrix_fp) from ch_determine_coding_consequences

  output:
  file "*consequences.tsv"

  script:
  """
  determine_coding_consequences.py --allele_fp ${allele_matrix_fp} --reference ${reference} --replicon ${replicon_id} > ${replicon_id}_alleles_var_cons0.95_consequences.tsv
  """
}


// Infer phylogeny
process infer_phylogeny {
  publishDir "${output_dir}", saveAs: { filename -> "${reference.simpleName}_${filename}" }

  input:
  tuple replicon_id, file(alignment_fp) from ch_infer_phylogeny

  output:
  file '*.tree'

  script:
  """
  FastTree -gtr -gamma -nt ${alignment_fp} > ${replicon_id}_alleles_var_cons0.95.tree
  """
}
*/
