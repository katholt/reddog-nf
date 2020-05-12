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


// Execute a command and capture standard streams along with return code
def execute_command(String command) {
  stdout = new StringBuilder()
  stderr = new StringBuilder()
  process = command.execute()
  process.waitForProcessOutput(stdout, stderr)
  return [process.exitValue(), stdout, stderr]
}


// For an optional stage param variable, check that it is either a Boolean or String
// If it is a string and either 'true' or 'false', return the boolean equivalent
def check_boolean_option(Object option, String name) {
  if (option.getClass() == java.lang.Boolean) {
    return option
  } else if (option.getClass() == java.lang.String) {
    if (option.toLowerCase() == 'true') {
      return true
    } else if (option.toLowerCase() == 'false') {
      return false
    }
  }
  exit 1, "error: ${name} option must be true or false"
}


// Check required input and outputs
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


// Require optional stage variables to be boolean
// We must check and change values if needed. The param variables are immutable so instead we declare new ones
run_read_subsample = check_boolean_option(params.subsample_reads, 'subsample_reads')
run_quality_assessment = check_boolean_option(params.quality_assessment, 'quality_assessment')
run_phylogeny = check_boolean_option(params.force_tree, 'force_tree')


// Check integer params
if (run_read_subsample & params.subsample_read_count.getClass() != java.lang.Integer) {
  if (! params.subsample_read_count.isInteger()) {
    exit 1, "error: subsample_read_count must be an integer, got {params.subsample_read_count}"
  }
}


// Create file objects from input parameters
reference_gbk_fp = file("${params.reference}")
output_dir = file("${params.output_dir}")


// Get name of reference
reference_name = reference_gbk_fp.simpleName


// Create channel for input read sets and get number of input isolates
Channel.fromFilePairs(params.reads, flat: true).into { ch_read_sets; ch_read_sets_count }
isolate_count = ch_read_sets_count.count()


// Ensure we have the appropriate input files
if (isolate_count.val < 1) {
  exit 1, "error: did not find any read files value '${params.reads}'"
}
if (! reference_gbk_fp.exists()) {
  exit 1, "error: reference input '${reference_gbk_fp}' does not exist"
}


// Set path to project bin directory
bin_dir = "${workflow.projectDir}/bin"


// Validate reference
script = 'validate_reference.py'
validate_command = "${bin_dir}/${script} --reference_fp ${params.reference}"
(return_code, stdout, stderr) = execute_command(validate_command)
if (return_code != 0) {
  exit 1, "error: validation of reference failed and exited with the following message:\n\n${stderr}"
}


// Validate first 10 readsets
script = 'validate_reads.py'
Channel.fromFilePairs(params.reads, flat: true).take(10).flatMap { it -> it[1..2] }.collect().set { ch_read_sets_validate }
validate_command = "${bin_dir}/${script} --reads_fps ${ch_read_sets_validate.val.join(' ')}"
(return_code, stdout, stderr) = execute_command(validate_command)
if (return_code != 0) {
  exit 1, "error: validation of first ten reads failed and exited with the following message:\n\n${stderr}"
}


// Do not run if output exists and contains files other than the run info directory (which is created by this point)
output_dir_files = []
output_dir.eachFile { output_dir_files.add(it.name) }
run_info_dirname = file(params.run_info_dir).simpleName
output_dir_files.remove(run_info_dirname)
if (output_dir_files.size() > 0 && ! params.force) {
  exit 1, "error: output directory '${output_dir}' already exists and contains other files, remove or use --force to overwrite"
}


// Do not run on MASSIVE unless user specifies profile to use to avoid inadvertently using a local executor
massive_hostnames = ['m3-login1', 'm3-login2']
on_massive = massive_hostnames.contains(InetAddress.getLocalHost().getHostName())
profile_explicit = workflow.commandLine.tokenize(' ').contains('-profile')
if (on_massive && ! profile_explicit) {
  exit 1, "error: to run on MASSIVE you must explicitly set -profile"
}


// Seed channels with data, allowing for optional subsampling
// Optional execution of subsampling is achieved by:
//   - populating only the subsample channel with input reads
//   - then feeding outputs of the subsample process to alignment and qc channels
//     - the `mix` operator is required to do this
//   - in the absence of subsampling, alignment and qc channels are populated as normal
ch_read_sets_qc = Channel.empty()
ch_read_sets_align = Channel.empty()
ch_read_sets_subsample = Channel.empty()
if (! run_read_subsample) {
  ch_read_sets.into { ch_read_sets_align; ch_read_sets_qc }
} else {
  ch_read_sets.set { ch_read_sets_subsample }
}


// Subsample reads if requested
process subsample_reads {
  input:
  tuple isolate_id, path(reads_fwd), path(reads_rev) from ch_read_sets_subsample

  output:
  tuple isolate_id, path('*1_subsampled*'), path('*2_subsampled*') into ch_subsampled_read_sets_align, ch_subsampled_read_sets_qc

  script:
  """
  subsample_reads.py --forward_fp ${reads_fwd} --reverse_fp ${reads_rev} --subsample_reads ${params.subsample_read_count} --output ./
  """
}


// Run read quality assessment
process create_read_quality_reports {
  publishDir "${output_dir}/fastqc/individual_reports/"

  input:
  tuple isolate_id, path(reads_fwd), path(reads_rev) from ch_read_sets_qc.mix(ch_subsampled_read_sets_qc)

  output:
  path '*html'
  path '*zip' into ch_fastqc_reports

  when:
  run_quality_assessment

  script:
  """
  fastqc -o . -t 1 -q -f fastq ${reads_fwd} ${reads_rev}
  """
}


// Aggregate reports with multiqc
process aggregate_read_quality_reports {
  publishDir "${output_dir}/fastqc/"

  input:
  path fastqc_report_fps from ch_fastqc_reports.collect()

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
  output:
  path '*.fasta' into ch_reference_fasta
  path '*fasta.*bt2' into ch_reference_bt2_index
  path '*fasta.fai' into ch_reference_samtools_index

  script:
  reference_fp = reference_name + '.fasta'
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
ch_reference_fasta.into { ch_align_reference; ch_call_snps_reference; ch_allele_matrix_reference }


// Align reads and apply some post-processing
//   - align with bowtie2
//     - minimum fragment length set to 2000, default of 500 too restrictive
//     - write metrics to output file, limit reporting to 30 minute intervals and final summary
//   - filtered unmapped reads
//   - sort reads
//   - index resulting BAM file
process align_reads_to_reference {
  publishDir "${output_dir}/bams/", pattern: '*.bam'

  input:
  tuple isolate_id, path(reads_fwd), path(reads_rev) from ch_read_sets_align.mix(ch_subsampled_read_sets_align)
  path reference_fp from ch_align_reference
  path reference_indices from ch_reference_bt2_index

  output:
  tuple val(isolate_id), path("${isolate_id}.bam"), path("${isolate_id}.bam.*") into ch_bams
  tuple val(isolate_id), path("${isolate_id}_metrics.tsv") into ch_alignment_metrics

  script:
  """
  bowtie2 --sensitive-local -x ${reference_fp} -1 ${reads_fwd} -2 ${reads_rev} -X 2000 --met 1800 --met-file ${isolate_id}_metrics.tsv | \
    samtools view -u -F4 | samtools sort > ${isolate_id}.bam
  samtools index ${isolate_id}.bam
  """
}
ch_bams.into { ch_call_snps_bams; ch_mpileup_bams; ch_mapping_stats_bams; ch_allele_matrix_bams }


// Create samtools multiway pileups
process create_mpileups {
  input:
  tuple isolate_id, path(bam_fp), path(bam_index_fp) from ch_mpileup_bams

  output:
  tuple val(isolate_id), path('*mpileup.tsv') into ch_mpileups
  path '*mpileup.tsv'  into ch_gene_coverage_depth_mpileups

  script:
  """
  samtools mpileup -aa ${bam_fp} > ${isolate_id}_mpileup.tsv
  """
}


// Get gene coverage and mean depth
process gene_coverage_depth {
  publishDir "${output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  path mpileup_fps from ch_gene_coverage_depth_mpileups.collect()

  output:
  path('gene_coverage.tsv')
  path('gene_depth.tsv')

  script:
  """
  create_coverage_depth_matrices.py --mpileup_fps ${mpileup_fps} --reference_fp ${reference_gbk_fp} --output_dir ./
  """
}


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
ch_call_snps_bams.join(ch_mpileups).set { ch_call_snps }

process call_snps {
  publishDir "${output_dir}/vcfs/", pattern: '*.vcf'

  input:
  tuple isolate_id, path(bam_fp), path(bam_index_fp), path(mpileup_fp) from ch_call_snps
  path reference_fp from ch_call_snps_reference
  path fasta_index from ch_reference_samtools_index

  output:
  tuple val(isolate_id), path('*_q30.vcf'), path('*_hets.vcf') into ch_snp_vcfs
  tuple val(isolate_id), path("${isolate_id}_coverage_depth.tsv") into ch_mapping_stats_coverage_depth
  path "${isolate_id}_sites.tsv" into ch_snp_sites_aggregate

  script:
  """
  # Call variants and index
  bcftools mpileup -Ou -f ${reference_fp} ${bam_fp} | bcftools call -Oz -c -v --ploidy 1 > ${isolate_id}_raw.bcf
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
ch_snp_vcfs.set { ch_mapping_stats_vcfs }


// Calculate mapping statistics
//   - get total reads from mapping metrics file
//   - coverage and mapped read counts previously calculated from mpileups
//   - get SNP, INDEL, and heterzygous SNP counts from VCFs
//   - calculate appropriate statisitics
//   - fail isolates on depth or coverage
//   - additionally for the largest replicon, require that at least 50% reads mapped
ch_mapping_stats_bams.join(ch_mapping_stats_vcfs)
  .join(ch_mapping_stats_coverage_depth)
  .join(ch_alignment_metrics)
  .set { ch_calculate_mapping_stats }

process calculate_mapping_statistics {
  input:
  tuple isolate_id, path(bam_fp), path(bam_index_fp), path(vcf_q30_fp), path(vcf_hets_fp), path(coverage_depth_fp), path(mapping_metrics_fp) from ch_calculate_mapping_stats

  output:
  path "${isolate_id}_mapping_stats.tsv" into ch_mapping_stats_aggregate

  script:
  """
  calculate_mapping_stats.py --bam_fp ${bam_fp} --vcf_q30_fp ${vcf_q30_fp} --vcf_hets_fp ${vcf_hets_fp} --coverage_depth_fp ${coverage_depth_fp} --mapping_metrics_fp ${mapping_metrics_fp} > ${isolate_id}_mapping_stats.tsv
  """
}


// Aggregate mapping statistics
//   - read in all records and aggregate by replicon
//   - for passing isolates, calculate phylogeny group using the ratio of 'total reads mapped' : 'coverage' / 100
//     - set maximum ratio for ingroup as 'ratio mean' + 'ratio stddev' * 2
//   - write aggregate stats to file, in the order of 'fail', 'outgroup', 'ingroup'
//   - get a list of isolates that pass mapping on at least one replicon
//   - get a list of replicons that have any passing isolate
//     - determined by having pass status and more than one SNP
process aggregate_mapping_statistics {
  publishDir "${output_dir}", pattern: '*_mapping_stats.tsv', saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  path(mapping_stats_fps) from ch_mapping_stats_aggregate.collect()

  output:
  path '*_mapping_stats.tsv'
  path 'isolate_replicons_passing.tsv' into ch_isolate_replicons_passing_filepath

  script:
  """
  aggregate_mapping_stats.py --rep_stats_fps ${mapping_stats_fps} --output_dir ./
  get_passing_isolate_replicons.awk ${mapping_stats_fps} > isolate_replicons_passing.tsv
  """
}


// Combine SNP sites
//   - create list of all unique high quality SNP sites
//   - exclude sites from isolates that failed for the respecitve replicon
process aggregate_snp_sites {
  input:
  path snp_sites_fps from ch_snp_sites_aggregate.collect()
  path isolate_replicons_passing_fp from ch_isolate_replicons_passing_filepath

  output:
  path 'snp_sites.tsv' into ch_snp_sites

  script:
  """
  aggregate_snp_sites.py --sites_fps ${snp_sites_fps} --isolate_replicons_passing_fp ${isolate_replicons_passing_fp} > snp_sites.tsv
  """
}


// Create allele matrix
//   - exclude isolates from that do not pass mapping to any replicon
//   - get genotype likelihoods for each SNP site with bcftool multiway pileup
//     - relaxed filtering to accept low-quality alleles at high quality sites
//   - for each replicon, write alleles called at each site as separate file

// Read in isolates the replicons that passed and generate a channel to emit [isolate_id, replicon_ids]
// Where replicon_ids is a string with each replicon_id separated by a single space
// Also get count of passing isolates so we can scale resource allocation if using SLURM executor
ch_isolate_replicons_passing_filepath.flatMap { filepath ->
    filepath.readLines().collect { line ->
      tokens = line.tokenize('\t')
      [tokens[0], tokens[1..-1].join(' ')]
    }
  }.into { ch_isolate_replicons_passing; ch_isolate_passing_count }
isolate_passing_count = ch_isolate_passing_count.count()

// Remove isolates that have no replicons that pass mapping criteria and then add SNP site file to each BAM
// We perform this here so that we do not run jobs for isolates that have no passing replicons
ch_allele_matrix_bams.join(ch_isolate_replicons_passing).combine(ch_snp_sites).set { ch_create_allele_matrix }

process create_allele_matrix {
  input:
  tuple isolate_id, path(bam_fp), path(index_fp), val(replicons_pass), path(sites_fp) from ch_create_allele_matrix
  path reference_fp from ch_allele_matrix_reference

  output:
  tuple val(isolate_id), path('*_alleles.tsv') into _ch_allele_matrix_aggregate

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
//   - filter merged matrix
//     - remove invariant sites; sites with only one allele (excluding '-')
//     - remove site if more than 5% of alleles are unknown (indicated by '-')

// Create a channel that emits allele matrices arranged by replicon_id
//   - input of [isolate_id, list(isolate_allele_matrices)]
//     - nextflow returns a list for multiple files or single object for one file
//     - check for different object types and process accordingly
//   - use isolate_id to robustly get replicon_id from allele matrix filename
//   - flat emit [replicon_id, isolate_allele_matrix] for each file
//   - group each matrix by replicon_id to emit [replicon_id, list(isolate_allele_matrices)]
_ch_allele_matrix_aggregate.flatMap { isolate_id, filepaths ->
    if (! (filepaths instanceof List)) {
      replicon_id = filepaths.getName().minus("_${isolate_id}_alleles.tsv")
      return [[replicon_id, filepaths]]
    } else {
      return filepaths.collect { filepath ->
          replicon_id = filepath.getName().minus("_${isolate_id}_alleles.tsv")
          [replicon_id, filepath]
        }
    }
  }.groupTuple().set { ch_allele_matrix_aggregate }

process aggregate_allele_matrices {
  publishDir "${output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  tuple replicon_id, path(allele_fps) from ch_allele_matrix_aggregate
  path snp_sites_fp from ch_snp_sites

  output:
  path '*_alleles.tsv'
  tuple val(replicon_id), path('*_alleles_core.tsv') into ch_allele_matrix

  script:
  """
  { head -n1 ${snp_sites_fp}; grep -w "^${replicon_id}" snp_sites.tsv; } > replicon_sites.tsv
  aggregate_allele_matrices.py --allele_fps ${allele_fps} --sites_fp replicon_sites.tsv > ${replicon_id}_alleles.tsv
  filter_allele_matrix.py --allele_fp ${replicon_id}_alleles.tsv > ${replicon_id}_alleles_core.tsv
  """
}
// Filter matrices that have no alleles so we don't needlessly execute downstream processes
ch_allele_matrix.filter { replicon_id, fp ->
    // Read first two lines of allele matrix and determine if we have data
    has_alleles = true
    fh = fp.newReader()
    for (int i = 0; i < 2; i++) { has_alleles = fh.readLine() != null }
    return has_alleles
  }.into { ch_snp_alignment_alleles; ch_consequences_alleles }


// Create SNP alignment
process create_snp_alignment{
  publishDir "${output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  tuple replicon_id, path(allele_matrix_fp) from ch_snp_alignment_alleles

  output:
  tuple val(replicon_id), path('*.mfasta') into ch_snp_alignment

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
  publishDir "${output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  tuple replicon_id, path(allele_matrix_fp) from ch_consequences_alleles

  output:
  path '*consequences_core.tsv'

  script:
  """
  determine_coding_consequences.py --allele_fp ${allele_matrix_fp} --reference_fp ${reference_gbk_fp} --replicon ${replicon_id} > ${replicon_id}_consequences_core.tsv
  """
}


// Infer phylogeny
process infer_phylogeny {
  publishDir "${output_dir}", saveAs: { filename -> "${reference_name}_${filename}" }

  input:
  tuple replicon_id, path(alignment_fp) from ch_snp_alignment

  output:
  path '*.tree'

  when:
  isolate_count.val <= 1000 | run_phylogeny


  script:
  """
  FastTree -gtr -gamma -nt ${alignment_fp} > ${replicon_id}_core.tree
  """
}
