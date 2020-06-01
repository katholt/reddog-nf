// Create allele matrix
//   - exclude isolates from that do not pass mapping to any replicon
//   - get genotype likelihoods for each SNP site with bcftool multiway pileup
//     - relaxed filtering to accept low-quality alleles at high quality sites
//   - for each replicon, write alleles called at each site as separate file
process create_allele_matrix {
  input:
  tuple isolate_id, path(bam_fp), path(index_fp), val(replicons_pass), path(sites_fp)
  path reference_fp

  output:
  tuple val(isolate_id), path('*_alleles.tsv'), emit: output

  script:
  """
  # Get sites only for replicons that pass mapping criteria
  { head -n1 ${sites_fp}; grep -wf <(echo ${replicons_pass} | tr ' ' '\n' | sed 's/^/^/') ${sites_fp}; } > isolate_sites.tsv
  # Create matrix
  create_allele_matrix.py --bam_fp ${bam_fp} --sites_fp isolate_sites.tsv --reference_fp ${reference_fp} --output_dir ./
  """
}


// Aggregate allele matrices
//   - group matrices by replicon id
//   - merge into single allele matrix
process aggregate_allele_matrices {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  tuple replicon_id, path(allele_fps)
  path snp_sites_fp
  val reference_name

  output:
  tuple val(replicon_id), path('*_alleles.tsv'), emit: output

  script:
  """
  { head -n1 ${snp_sites_fp}; grep -w "^${replicon_id}" snp_sites.tsv; } > replicon_sites.tsv
  aggregate_allele_matrices.py --allele_fps ${allele_fps} --sites_fp replicon_sites.tsv > ${replicon_id}_alleles.tsv
  """
}


// Filter allele matrices
//   - remove invariant sites; sites with only one allele (excluding '-')
//   - remove site if more than 5% of alleles are unknown (indicated by '-')
process filter_allele_matrix {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  tuple replicon_id, path(allele_fp)
  val reference_name

  output:
  tuple val(replicon_id), path('*_alleles_core.tsv'), emit: output

  script:
  """
  filter_allele_matrix.py --allele_fp ${allele_fp} --conservation ${params.allele_matrix_cons} > ${replicon_id}_alleles_core.tsv
  """
}


