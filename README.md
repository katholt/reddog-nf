# reddog-nf


## Processing differences
* Most stages have been rearranged and combined in some way
    - overall pipeline structure is preserved
* Validation on first few readsets and input reference
* Quality assessment of input reads
    - fastqc and aggregated by multiqc
* Reads subsampling if requested
    - useful to ensure read alignment can execute on shortq
* Isolates in replicon statistics file are ordered as failed, outgroup, ingroup
* The coding consequences process additionally provides a list of affected isolates
* Filtering BAMs in single-shot while mapping
    - unmapped read statistic taken from bowtie2 metrics file
    - no unfiltered BAM intermediate file
* Using `bcftools mpileup` rather than `samtools mpileup` for variant calling
    - `bcftools mpileup` has built in filtering criteria causing difference in numbers of SNP calls
* Using Kat's approach to determine alleles at high quality SNP sites
    - this impacts allele filtering at 95% conservation as alleles can differ from previous method
    - usually less sites with known alleles as we now only consider >=Q20 alignments for consensus calls
* During matrix aggregation ingroup/outgroup is not calculated for replicons with < 2 passing isolates
    - instead the group status is set to undetermined
* Coding consequences does not assess alleles in features that have compound locations
    - occurs when gene crosses contig boundary or due to annotation
    - several genes in the Hi reference have compound location annotations to just show frameshifts
        - determining consequences here is not feasible without defining what is the true ORF
* No equivalent for following files
    - gene presence/absence
    - gene counts in isolate set
    - AllStats


## TODO
* Remove CheckInput and CheckOutput classes
* Further develop read simulator for testing
    - currently simulates SNPs, INDELs, and low quality alleles
        - where low quality alleles appear as '-' in the allele table
    - junk reads
        - tests for mapping stats and ingroup and outgroup assignment
        - best approach would be to adjust coverage to make room for junk reads
            - conditionally reduce step appropriately to get lower good read count
            - add junk reads to achieve desired coverage/ read count
        - best approach will required some refactor to keep code clean
            - classes are probably needed at this point tbh
* Automate comparison of test output data
* FastQC seems to give wrong phred score for second isolate in read simulation run
* I think read validation fails if format is correct but there is very short reads
    - probably want to abort if reads are so short mapping is not useful
    - accidentally gave nf reddog simulated reads with no sequence and it tried to run
* Subsampling approach is too slow when subsampling down to a large number
    - subsampling to 5 million reads takes ~10 minutes
    - likely a result of hash collisons during set operations
    - different approach could be to use a probability function at each read iteration
        - output read counts will not be exact but will scale much better
        - may or may not be quicker at lower numbers
* Model memory usage and wall time ~ readsets file size
    - the is a strongly relationship (but other factors are important)
    - will generally provide very good of run time in read alignment, SNP calling
    - use data from the nextflow run reports - very easy to extract
    - also look at resource usage for jobs that aggregate
        - model on Q30 SNPs, passing isolates, replicon size, etc
* Coding consequences currently only examines each SNP individually
    - report coding result when multiple SNPs are in one codon
    - this probably should replace looking at SNPs individually
    - would require reworking a large amount of `determine_coding_consequences.py`
* Scaling job resources on file size may be more appropriate
    - use the `size()` method, returns file size in bytes
* Set output files to be copied or moved rather than symlinked
    - move is preferable, particularly for large runs producing many BAMs
    - could do this through nextflow publishDir interface
        - if moving files this will cause problems unless they're not used downstream
    - or through a final stage to run some BASH to move files
        - for file in $(find . -type l); do mv $(readlink ${file}) ${file}; done
* Coding consequences for hets
    - create matrix of hets in the same way we do for homs
    - this should be done conditionally
    - pass the data forward
* Mixed sample detection
    - SNP:het ratio in replicon statistics
    - use zoe's scripts to plot
        - reads het and q30 vcfs
        - plot %mapped to alt, coloured by hets or homs
* Add option to send email of nextflow fail (pipeline failure, not task failure)
* Probably can remove mapped flag check in `get_reads_mapped.awk`
    - was previously passing unfiltered bam
* For `get_snp_sites.awk`, check that the input ref and alt allele can be the same (i.e. not use of '.')
* Read simulator in `bin/` and be removed for distribution


## Planned improvements
* Run report with additional information, see RedDog and Jane's pipeline
* Use recommended filtering for variant filtering in samtools
    - https://samtools.github.io/bcftools/howtos/filtering.html
* Currently for fail samples, the largest replicon requires 50% of reads mapped
    - we should do a proportional requirement i.e. require n% ~ replicon\_size
    - simple additional to `calculate_replicon_statistics.py`
* Improve allele matrix creation - for each isolate we're writing both position and reference
    - at least half the disk i/o by removing reference column in each
* Chunked reads for SNP aligments and matrix aggregation
    - avoid consuming very large amounts of memory for large datasets
* Allow FASTA input
    - must disable certain processes - use `when` directive if clean enough
    - we'll probably be able to use a single `when` directive early in the pipeline
* Sort inputs by file size
    - can provide small reduction in runtime


## Planned significant additions
* Merge run pipeline
    - separate nextflow script
* Single reads pipeline
    - separate nextflow script
* Variants other than SNPs


## Queries
* Is a SNP site defined as variant if all isolates have the same allele that is different to reference
    - presumably yes but should check
* During variant calling there are two quality filters
    - vcfutils varFilter: >= 30 RMS mapping quaility
    - python script: >= 30 mapping quality
    - these are different statistics but is the raw mapping quality filter sufficient?
* Do we need unmapped reads?
    - currently removing all unmapped reads from BAMs
    - this could be produced optionally
* Counts of genes in isolate set with >=n depth/coverage
    - this would be done in the gene\_coverage\_depth process
    - trival to calculate outside of pipeline
    - check if this is something wanted
* Ingroup defined as:
    - m = total_reads / replicon_coverage / 100
    - isolate's m <= (stddev all isolate's m) * 2
    - should we also look at number of SNPs?
* Hets are currently determined by looking at a diploid genotype bcftools call
    - specifically looking at GT fields of 1/1 or 0/0
    - what are the implications of callign in diploid mode rather than haploid?
    - we can look at variant type and AF1 to achieve the same in haploid, as far as I can tell
    - also SNPs with ONE other non-reference allele are considered hets
        - e.g. ref: A; alt: T (50 reads), C (1 read)
        - should this be a hom of A->T rather than a het?
* Other mixed sample detection methods?


## Misc notes
* I observed the pipeline to continue after exceeding maximum retries for a job
    - occured on haemophilus dataset
    - was not able to recreate on small test pipeline


## Items for further testing
This is a list of processes/outputs that are yet to be closely tested
* New allele matrix creation method
    - are call sites correct?
* Allele matrix filtering, large dataset with many unknown SNPs
    - exclusion of SNP positions which lack sufficient data
* Ingroup/output calling, large dataset with many unknown SNPs
    - ratio calculation and equality comparison
* Output of gene\_coverage\_depth process
* SNP alignment
* Phylogeny
* Coding consequence
    - interval tree balance
    - ensure bounds capture
    - consider comparing against https://samtools.github.io/bcftools/howtos/csq-calling.html if applicable
* Isolate list in consequences containing more than one isolate
* Coding consequences for genes that join across contig boundaries (on same contig)
