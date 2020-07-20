# Important note
* Not currently using a strand bias cutoff
    - this is only used in reddog when the runType is set to PE
    - though runType is only ever set to phylogeny or pangenome
    - it may be that runType should be replaced with readType in rr

# TODO
* Currently copying merge data outside of nextflow execution
    - for very large datasets, the script may appear to hang during this operation
    - we could notify user of the operation w progress?
    - or perhaps create a process just to perform the copy w the 'local' executor
* Rewrite gene coverage/depth process to be more efficient
    - see report.html for flexneri
* Investigate use of async non-nf process code so we can:
    * Validate some input readsets
    * Warn when read quality assessment is requested on merge run but merge has fastqc data
* Coding consequences for hets
    - create matrix of hets in the same way we do for homs
    - this should be done conditionally
    - pass the data forward
* Mixed sample detection
    - SNP:het ratio in replicon statistics
    - use zoe's scripts to plot
        - reads het and q30 vcfs
        - plot %mapped to alt, coloured by hets or homs
* Investigate scope access of processes for purpose of variable passing
    - I would also like to use variables or value channels set within workflows in process closures
        - e.g. for process resources: time = { (isolates_passing.val * 0.5).MB * task.attempt }
        - this was possible in DSL1 but I'm yet to find an equivalent in DSL2
    - rn we're passing reference name by argument very often
    - this is a little messy tbh, hopefuly can do better
* MultiQC/FastQC seems to give wrong phred score for second isolate in read simulation run


## Queries
* During variant calling there are two quality filters
    - vcfutils varFilter: >= 30 RMS mapping quaility
    - python script: >= 30 mapping quality
    - these are different statistics but is the raw mapping quality filter sufficient?
* Ingroup defined as:
    - m = total\_reads / replicon\_coverage / 100
    - isolate's m <= (stddev all isolate's m) * 2
    - should we also look at number of SNPs?
* Hets are currently determined by looking at a diploid genotype bcftools call
    - specifically looking at GT fields of 1/1 or 0/0
    - what are the implications of callign in diploid mode rather than haploid?
    - we can look at variant type and AF1 to achieve the same in haploid, as far as I can tell
    - also SNPs with ONE other non-reference allele are considered hets
        - e.g. ref: A; alt: T (50 reads), C (1 read)
        - should this be a hom of A->T rather than a het?


## Items to be closely tested
* Output of gene\_coverage\_depth process
* SNP alignment and phylogeny
