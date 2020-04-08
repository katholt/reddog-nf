# reddog-nf


## Processing differences
* No checkBam stage equivalent yet
    - this will be a simple post-process script
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


## Currently known output differences
* New versions of software calling additional SNPs with same parameters
    - Presumably this is due to improvements and bug fixes


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
