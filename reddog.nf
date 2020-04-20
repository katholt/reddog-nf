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
reference_gbk_fp = file("${params.reference}")
output_dir = file("${params.output_dir}")


// Get name of reference
reference_name = reference_gbk_fp.simpleName


// Create channel for input read sets
Channel.fromFilePairs(params.reads, flat:true).set { ch_read_sets }


// Prepare reference
//   - convert to FASTA
//   - index FASTA with bowtie2 and samtools
//   - get list of replicon names and check they're unique
process prepare_reference {
  output:
  file '*.fasta' into ch_reference_fasta
  file '*fasta.*bt2' into ch_reference_bt2_index
  file '*fasta.fai' into ch_reference_samtools_index

  script:
  reference_fp = reference_gbk_fp.simpleName + '.fasta'
  """
  # Convert to FASTA format
  genbank_to_fasta.py --input_fp ${reference_gbk_fp} --output_fp ${reference_fp}

  # Get replicon names then check they are all unique
  replicons=\$(grep '^>' ${reference_fp} | sed -e 's/^>//' -e 's/ .\\+\$//')
  replicons_duplicate=\$(echo "\${replicons}" | sort | uniq -d)
  replicons_duplicate_count=\$(echo -n "\${replicons_duplicate}" | sort | uniq -d | wc -l)
  if [[ "\${replicons_duplicate_count}" -gt 0 ]]; then
    echo "Found duplicate replicon names, you must rename these" 1>&2
    echo "Duplicates: \${replicons_duplicate}" 1>&2
    exit 1
  fi

  # Index
  bowtie2-build ${reference_fp} ${reference_fp}
  samtools faidx ${reference_fp}
  """
}
ch_reference_fasta.into {
    ch_align_reference;
    ch_call_snps_reference;
    ch_allele_matrix_reference;
  }


// Align reads and apply some post-processing
//   - align with bowtie2
//     - minimum fragment length set to 2000, default of 500 too restrictive
//     - write metrics to output file, limit reporting to 30 minute intervals and final summary
//   - filtered unmapped reads
//   - sort reads
//   - index resulting BAM file
process align_reads_to_reference {
  input:
  tuple isolate_id, file(reads_fwd), file(reads_rev) from ch_read_sets
  file reference_fp from ch_align_reference
  file reference_indices from ch_reference_bt2_index

  output:
  tuple val(isolate_id), file("${isolate_id}.bam"), file("${isolate_id}.bam.*") into ch_bams
  tuple val(isolate_id), file("${isolate_id}_metrics.tsv") into ch_alignment_metrics

  script:
  """
  bowtie2 --sensitive-local -x ${reference_fp} -1 ${reads_fwd} -2 ${reads_rev} -X 2000 --met 1800 --met-file ${isolate_id}_metrics.tsv | \
    samtools view -u -F4 | samtools sort > ${isolate_id}.bam
  samtools index ${isolate_id}.bam
  """
}
ch_bams.into { ch_call_snps_bams; ch_replicon_stats_bams; ch_allele_matrix_bams }


// Call SNPs and record high quality SNP sites
//   - get genotype likelihoods with a multiway pileup
//   - call variants in consensus caller mode
//   - filter variants - additional to defaults, require minimum depth of 5 and minimum RMS mapping quality of 30
//   - split variants into homozygous SNPs and heterzygous SNPs
//     - require all variants to have a mapping quality of 30 (different to RMS mapping quality)
//     - homozygous if allele frequency of 1 or a genotype different to reference, else heterozygous if not an INDEL
//   - get high quality SNP sites
process call_snps {
  input:
  tuple isolate_id, file(bam_fp), file(bam_index_fp) from ch_call_snps_bams
  file reference_fp from ch_call_snps_reference
  file fasta_index from ch_reference_samtools_index

  output:
  tuple val(isolate_id), file('*_q30.vcf'), file('*_hets.vcf') into ch_snp_vcfs
  file "${isolate_id}_sites.tsv" into ch_snp_sites_aggregate

  script:
  """
  bcftools mpileup -Ou -f ${reference_fp} ${bam_fp} | bcftools call -cv --ploidy 1 | \
    vcfutils.pl varFilter -d5 -Q30 | filter_variants.awk -v snp_fp=${isolate_id}_q30.vcf -v het_fp=${isolate_id}_hets.vcf
  get_snp_sites.awk ${isolate_id}_q30.vcf > ${isolate_id}_sites.tsv
  """
}
ch_snp_vcfs.set { ch_replicon_stats_vcfs }


// Calculate replicon statistics
//   - get total reads from mapping metrics file
//   - get coverage and mapped read counts from samtools mpileup and samtools view
//   - get SNP, INDEL, and heterzygous SNP counts from VCFs
//   - calculate appropriate statisitics
//   - fail isolates on depth or coverage
//   - additionally for the largest replicon, require that at least 50% reads mapped
ch_replicon_stats_bams.join(ch_replicon_stats_vcfs).join(ch_alignment_metrics).set { ch_calculate_replicon_stats }
process calculate_replicon_statistics {
  input:
  tuple isolate_id, file(bam_fp), file(bam_index_fp), file(vcf_q30_fp), file(vcf_hets_fp), file(mapping_metrics_fp) from ch_calculate_replicon_stats

  output:
  file "${isolate_id}_replicon_stats.tsv" into ch_replicon_stats_aggregate

  script:
  """
  calculate_replicon_stats.py --bam_fp ${bam_fp} --vcf_q30_fp ${vcf_q30_fp} --vcf_hets_fp ${vcf_hets_fp} --mapping_metrics_fp ${mapping_metrics_fp} > ${isolate_id}_replicon_stats.tsv
  """
}


// Aggregate replicon statistics
//   - read in all records and aggregate by replicon
//   - for passing isolates, calculate phylogeny group using the ratio of 'total reads mapped' : 'coverage' / 100
//     - set maximum ratio for ingroup as 'ratio mean' + 'ratio stddev' * 2
//   - write aggregate stats to file, in the order of 'fail', 'outgroup', 'ingroup'
//   - get a list of samples that pass mapping on at least one replicon
//   - get a list of replicons that have any passing isolate
//     - determined by having pass status and more than one SNP
process aggregate_replicon_statistics {
  input:
  file(replicon_stats_fps) from ch_replicon_stats_aggregate.collect()

  output:
  file "*_RepStats.tsv" into ch_replicon_stats
  env samples_pass into _ch_samples_pass
  env replicons_pass into _ch_replicons_pass

  script:
  """
  aggregate_replicon_stats.py --rep_stats_fps ${replicon_stats_fps} --output_dir ./
  samples_pass=\$(grep -wh 'p' *_RepStats.tsv | cut -f1 -d\$'\t' | sort | uniq)
  replicons_pass=\$(awk '\$10 == "p" && \$7 > 0 { print FILENAME }' *RepStats.tsv | sort | uniq | sed 's/_RepStats.tsv//')
  """
}
ch_replicon_stats.into { ch_snp_sites_replicon_stats; ch_allele_matrix_replicon_stats }
_ch_samples_pass.tokenize(' ').flatMap().set { ch_samples_pass }
_ch_replicons_pass.tokenize(' ').flatMap().set { ch_replicons_pass }


// Combine SNP sites
//   - create list of all unique high quality SNP sites
//   - exclude sites from isolates that failed for the respecitve replicon
process get_snp_sites {
  input:
  file snp_sites_fps from ch_snp_sites_aggregate.collect()
  file replicon_stats_fps from ch_snp_sites_replicon_stats

  output:
  file "snp_sites.tsv" into ch_snp_sites

  script:
  """
  aggregate_snp_sites.py --sites_fps ${snp_sites_fps} --replicon_stats_fp ${replicon_stats_fps} > snp_sites.tsv
  """
}


// Create allele matrix
//   - exclude isolates from that do not pass mapping to any replicon
//   - get genotype likelihoods for each SNP site with bcftool multiway pileup
//     - relaxed filtering to accept low-quality alleles at high quality sites
//   - for each replicon, write alleles called at each site as separate file

// Remove isolates that have no replicons that pass mapping criteria and then add SNP site file to each BAM
ch_samples_pass.join(ch_allele_matrix_bams).combine(ch_snp_sites).set { ch_create_allele_matrix }

process create_allele_matrix {
  input:
  tuple isolate_id, file(bam_fp), file(index_fp), file(sites_fp) from ch_create_allele_matrix
  file replicon_stats_fps from ch_allele_matrix_replicon_stats
  file reference_fp from ch_allele_matrix_reference

  output:
  file '*_alleles.tsv'  into _ch_allele_matrix_aggregate

  script:
  """
  # Get sites only for replicons that pass mapping criteria
  replicons_pass=\$(grep -l NCTC13753_set2 *RepStats.tsv | sed 's/_RepStats.tsv//')
  { head -n1 ${sites_fp}; grep -wf <(echo \${replicons_pass} | tr ' ' '\n' | sed 's/^/^/') ${sites_fp}; } > isolate_sites.tsv
  # Create matrix
  create_allele_matrix.py --bam_fp ${bam_fp} --sites_fp isolate_sites.tsv --reference_fp ${reference_fp} --output_dir ./
  """
}


// Aggregate allele matrices
_ch_allele_matrix_aggregate.collect().toList().combine(ch_replicons_pass).set { ch_allele_matrix_aggregate }
process aggregate_allele_matrices {
  input:
  tuple replicon_id, file(allele_fps) from ch_allele_matrix_aggregate

  output:
  file('*alleles_var*tsv')

  script:
  """
  # Aggregate matrices then filter
  #filter_allele_matrix.py --allele_fp aggregate_alleles.tsv > ${replicon_id}_alleles_var_cons0.95.tsv
  echo ${replicon_id}
  echo ${allele_fps}
  touch temp_alleles_var.tsv
  exit 1
  """
}


/*
// Create SNP alignment
process create_snp_alignment{
  publishDir "${output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

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
  publishDir "${output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

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
  publishDir "${output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

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
