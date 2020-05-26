// Call variants, filter for SNPs, and record high quality SNP positions
//   - get genotype likelihoods with a multiway pileup
//   - call variants in consensus caller mode
//   - get average depth of each replicon
//   - filter variants - additional to defaults, require:
//     - minimum depth of 5
//     - minimum RMS mapping quality of 30
//     - maximum depth of mean replicon depth * 2
//   - split variants into homozygous SNPs and heterzygous SNPs
//     - require all variants to have a mapping quality of 30 (different to RMS mapping quality)
//     - homozygous if allele frequency of 1 or a genotype different to reference, else heterozygous if not an INDEL
//   - get high quality SNP sites
process call_snps {
  publishDir "${params.output_dir}/vcfs/", pattern: '*.vcf'

  input:
  tuple isolate_id, path(bam_fp), path(bam_index_fp), path(mpileup_fp)
  path reference_fp
  path fasta_index

  output:
  tuple val(isolate_id), path('*_q30.vcf'), path('*_hets.vcf'), emit: vcfs
  tuple val(isolate_id), path("${isolate_id}_coverage_depth.tsv"), emit: coverage_depth
  path "${isolate_id}_sites.tsv", emit: sites

  script:
  """
  # Call variants and index
  bcftools mpileup -Ou -f ${reference_fp} ${bam_fp} | bcftools call -Oz -c -v > ${isolate_id}_raw.bcf
  bcftools index ${isolate_id}_raw.bcf
  # Get depths for filtering purposes
  get_coverage_depth.awk ${mpileup_fp} > ${isolate_id}_coverage_depth.tsv
  # Filter SNPs for each replicon with appropriate max depth
  while read -a data; do
    replicon_id=\${data[0]};
    depth_max=\${data[1]};
    bcftools view -r \${replicon_id} ${isolate_id}_raw.bcf | vcfutils.pl varFilter -d5 -D\${depth_max} -Q30 > ${isolate_id}_\${replicon_id}.vcf
  done < <(awk '{print \$1, \$2 * 2}' ${isolate_id}_coverage_depth.tsv)
  cat ${isolate_id}*vcf | filter_variants.awk -v snp_fp=${isolate_id}_q30.vcf -v het_fp=${isolate_id}_hets.vcf
  # Get SNP sites
  get_snp_sites.awk ${isolate_id}_q30.vcf > ${isolate_id}_sites.tsv
  """
}
