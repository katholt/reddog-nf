// Subsample se reads
process subsample_reads_se {
  input:
  tuple isolate_id, path(reads_fp)

  output:
  tuple isolate_id, path('*subsampled.fastq.gz'), emit: output

  script:
  """
  seqtk sample -s0 ${reads_fp} ${params.subsample_read_count} | gzip > ${isolate_id}_subsampled.fastq.gz
  """
}


// Subsample pe reads
process subsample_reads_pe {
  input:
  tuple isolate_id, path(reads_fwd), path(reads_rev)

  output:
  tuple isolate_id, path('*_subsampled_1.fastq.gz'), path('*_subsampled_2.fastq.gz'), emit: output

  script:
  """
  seqtk sample -s0 ${reads_fwd} ${params.subsample_read_count} | gzip > ${isolate_id}_subsampled_1.fastq.gz
  seqtk sample -s0 ${reads_rev} ${params.subsample_read_count} | gzip > ${isolate_id}_subsampled_2.fastq.gz
  """
}


// Run read quality assessment
process create_read_quality_report {
  publishDir "${params.output_dir}/fastqc/individual_reports/", mode: 'copy'

  input:
  path reads_fp

  output:
  path '*html'
  path '*zip', emit: output

  script:
  """
  fastqc -o . -t 1 -q -f fastq ${reads_fp}
  """
}


// Aggregate reports with multiqc
process aggregate_read_quality_reports {
  publishDir "${params.output_dir}/fastqc/", mode: 'copy'

  input:
  path fastqc_report_fps

  output:
  path '*html'
  path 'multiqc_data/'

  script:
  graph_type_opt = fastqc_report_fps.size() < 100 ? '--interactive' : '--flat'
  """
  multiqc . -m fastqc --data-dir ${graph_type_opt}
  """
}

// Prepare reference
//   - convert to FASTA
//   - index FASTA with bowtie2 and samtools
//   - get list of replicon names and check they're unique
process prepare_reference {
  executor 'local'

  input:
  file reference_gbk_fp

  output:
  val reference_gbk_fp.simpleName, emit: name
  path '*.fasta', emit: fasta
  path '*fasta.*bt2', emit: bt2_index
  path '*fasta.fai', emit: samtools_index

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


// Create samtools multiway pileups
process create_mpileups {
  input:
  tuple isolate_id, path(bam_fp), path(bam_index_fp)

  output:
  tuple val(isolate_id), path('*mpileup.tsv'), emit: output
  path '*mpileup.tsv', emit: output_noid

  script:
  """
  samtools mpileup -aa ${bam_fp} > ${isolate_id}_mpileup.tsv
  """
}


// Combine SNP sites
//   - create list of all unique high quality SNP sites
//   - exclude sites from isolates that failed for the respecitve replicon
process aggregate_snp_sites {
  input:
  path snp_sites_fps
  path replicons_passing_fp

  output:
  path 'snp_sites.tsv', emit: output

  script:
  """
  aggregate_snp_sites.py --sites_fps ${snp_sites_fps} --replicons_passing_fp ${replicons_passing_fp} > snp_sites.tsv
  """
}


// Create SNP alignment
process create_snp_alignment {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }, mode: 'copy'

  input:
  tuple replicon_id, path(allele_matrix_fp)
  val reference_name

  output:
  tuple val(replicon_id), path('*.mfasta'), emit: output

  script:
  """
  create_snp_alignment.py --allele_fp ${allele_matrix_fp} > ${replicon_id}_core.mfasta
  """
}


// Determine coding consequences
//   - create interval tree from gene features
//   - for each site find genes that it falls within
//   - infer codon change and amino acid change
process determine_coding_consequences {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }, mode: 'copy'

  input:
  tuple replicon_id, path(allele_matrix_fp)
  path reference_gbk_fp

  output:
  path '*consequences_core.tsv'

  script:
  reference_name = reference_gbk_fp.simpleName
  """
  determine_coding_consequences.py --allele_fp ${allele_matrix_fp} --reference_fp ${reference_gbk_fp} --replicon ${replicon_id} > ${replicon_id}_consequences_core.tsv
  """
}


// Infer phylogeny
process infer_phylogeny {
  publishDir "${params.output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }, mode: 'copy'

  input:
  tuple replicon_id, path(alignment_fp)
  val reference_name
  val pass_count
  val run_phylogeny

  output:
  path '*.tree'

  when:
  pass_count <= 1000 | run_phylogeny

  script:
  """
  FastTree -gtr -gamma -nt ${alignment_fp} > ${replicon_id}_core.tree
  """
}
