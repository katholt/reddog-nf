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


## TODO
* Provide mean replicon mapping depth (or threshold) to variant calling process
    - average depth * 2 for variant filtering
    - some inline groovy code will be required to assign values to appropriate replicons
* Coding consequences for hets
    - this will take some work
    - allele information needs to be obtained, represent as sparse matrix?
      - or as full matrix but with 'placeholder' values for homozygous SNPs
    - new process/ python script to run code - must explicity indicate which isolate has the het
      - a het in one isolate doesn't imply the site is heterozygous in all isolates
* Conditionally execute phylogeny process
    - on the basis of having n or more isolates
* Check we have input read sets
    - create different channel to check the first item
* Currently for fail samples, the largest replicon requires 50% of reads mapped
    - we should do a proportional requirement i.e. require n% ~ replicon_size
    - simple additional to `calculate_replicon_statistics.py`
* Improve allele matrix creation - for each isolate we're writing both position and reference
    - at least half the disk i/o by removing reference column in each
* Chunked reads for SNP aligments and matrix aggregation
    - avoid consuming very large amounts of memory for large datasets
* Replace excessive use of list, dict, set comprehensions with more digestible for loops
    - see `calculate_replicon_statisitics.py` for the main offender
* Probably can remove mapped flag check in get_reads_mapped.awk
    - was previously passing unfiltered bam
* For get_snp_sites.awk, check that the input ref and alt allele can be the same (i.e. not use of '.')
* Revert min. base quality to 20 in the site-specific consensus calling of create_allele_matrix.py
    - reduced required quality to use with small test dataset
* Run report with additional information, see RedDog and Jane's pipeline
* Add optional fastqc
    - default behaviour to generate reports
    - option to turn off
* Mixed sample detection
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
