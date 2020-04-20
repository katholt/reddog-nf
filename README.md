# reddog-nf


## Processing differences
* The calculate_replicon_statistics process
    - combines several steps
        - getSamStats
        - getCoverage
        - deriveRepStats
    - calculates the following from combination of samtools view/mpileup and awk scripts:
        - coverage and depth
        - heterozygous SNP count
        - SNPs and INDELs counts
        - read counts
* Isolates in replicon statistics file are ordered as failed, outgroup, ingroup
* The coding consequences process additionally provides a list of affected isolates
    - implemented an interval tree to quickly collect features containing a given SNP
* Filtering BAMs in single-shot while mapping
    - unmapped read statistic taken from bowtie2 metrics file
* Using `bcftools mpileup` rather than `samtools mpileup` for variant calling
* Not using maximum depth for vcfutils.pl varFilter
    - Reddog uses the 2 * average depth of each replicon
    - greatly simplifies processes and provides speed improvement, will work to reintroduce
* Coding consequences for hets (TODO)


## TODO
* Provide passing isolates + replicons to aggregate_allele_matrices
    - solution is to allow the python script to robustly glob the isolate allele matrix for a given replicon
    - give indirectly by replicon stats (and calculate) or directly by nf channel (this seems difficult rn)

* Check we have input read sets
    - create different channel to check the first item
* Likely need to add maximum depth for variant calls
    - use 2 * average depth
* Currently for fail samples, the largest replicon requires 50% of reads mapped
    - we should do a proportional requirement i.e. require n% ~ replicon_size
    - simple additional to `calculate_replicon_statistics.py`
* Probably can remove mapped flag check in get_reads_mapped.awk
    - was previously passing unfiltered bam
* Replace pseudo-genome alignment with SNP alignment
* Add gene coverage and presence/absence matrix
    - needs depth and coverage for each replicon
    - may be best to reintroduce getCoverage stage
        - output will be used by both this and the relicon stats process
* Remove processing per replicon
    - just do all at once and handle in scripts
    - additional code for replicon statistics, SNP sites, allele matrix
* Conditionally execute phylogeny process
    - on the basis of having n or more isolates
* For get_snp_sites.awk, check that the input ref and alt allele can be the same (i.e. not use of '.')
* Revert min. base quality to 20 in the site-specific consensus calling of create_allele_matrix.py
    - reduced required quality to use with small test dataset
* Parallelise process that are working on sample basis
    - achieve this with the `buffer(size: n)` operator and `parallel` within the process
    - this woud require unpacking of list items in the script block, a little messy
    - still not a complete solution, this will take some consideration
    - collect_snp_sites, create_allele_matrix, collect_snp_sites
* Parallelise sort in aggregate_snp_sites process
* Harden creation of the ch_sample_status channel
    - do some BASH comparisons, assertions to ensure that we have good data
* For create_allele_matrix, there may be some condition where no output files are created, check carefully
    - when a replicon has not hihg quality sites, it is skipped
    - when a sample has no high quality consensus calls, it will be skipped
    - however bcftools *may* always produce entries for sites even if a call cannot be made, check
* Remove dependency of consistent ordering between nf and BASH
    - done in transform of _ch_replicon_stats_vcfs for calculate_replicon_statistics
    - just get replicon_id from filename, this is standard practice
    - maybe I can later find a more optimal way later on
* I've used `Channel.toList().toList()` inorder to get desire behaviour from `combine`
    - is there a better way?
* Run report with additional information, see RedDog and Jane's pipeline
* Add optional fastqc
    - default behaviour to generate reports
    - option to turn off
* Mixed sampel detection
    - SNP:het ratio in replicon statistics?
    - use zoe's scripts to plot
        - reads het and q30 vcfs
        - plot %mapped to alt, coloured by hets or homs
* Allow FASTA input
    - must disable certain processes - use `when` directive if clean enough
    - we'll probably be able to use a single `when` directive early in the pipeline
* Check user inputs and configuration before executing pipeline
* Sort inputs by file size
    - can provide small reduction in runtime
* Mixed sample detection
    - to be discussed with the lab
* Merge run pipeline
    - separate nextflow script
* Single reads pipeline
    - separate nextflow script


## Items for further testing
This is a list of processes/outputs that are yet to be closely tested
* Allele matrix filtering, large dataset with many unknown SNPs
    - exclusion of SNP positions which lack sufficient data
* Ingroup/output calling, large dataset with many unknown SNPs
    - ratio calculation and equality comparison
* Pseudo-genome alignment
* Phylogeny
* Coding consequence
    - interval tree balance
    - ensure bounds capture
