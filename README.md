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
* Allele matrix creation is different, currently requires improvement
    - creates final matrix for each replicon in one script with no parallelisation
    - most important impact of this is the generation of the SNP position set
        - this is currently not tractable for large datasets, will be reworked next week
        - small batched jobs to create position set and individual matrices as in RedDog
* The coding consequences process additionally provides a list of affected isolates
    - implemented an interval tree to quickly collect features containing a given SNP


## TODO
* Fix filter_snps_vcfutils process
    - incorrect parameterisation of -D
* Update allele matrix creation
    - using Kat script
        - uses q30 variants for SNP positions
        - replaces all site consensus with q20 site-specific consensus
    - split into smaller jobs
        - batch into jobs of ~50-100 isolates
    - add aggregation step
        - specify current directory as input, select files with glob generator in python
* Replace pseudo-genome alignment with SNP alignment
* Add gene coverage and presence/absence matrix
    - needs depth and coverage for each replicon
    - may be best to reintroduce getCoverage stage
        - output will be used by both this and the relicon stats process
* Conditionally execute phylogeny process
    - on the basis of having n or more isolates
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
