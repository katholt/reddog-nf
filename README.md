# reddog-nf


## Processing differences
* Most stages have been rearranged and combined in some way
    - overall pipeline structure is preserved
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
* No equivalent for following files
    - gene presence/absence
    - gene counts in isolate set
    - AllStats


## TODO
* Pipeline continues when a job fails n times and generates output files
    - there is not immediate notification of this occuring
    - seems to be unwanted behaviour
    - dynamically setting errorStrategy seems to be a possible option:
        - errorStrategy: { task.attempt < 3 ? 'retry' : 'terminate' }
        - on third set errorStrategy to terminate
        - may be able to access the maxTries directive to do:
            - errorStrategy: { task.attempt < task.maxTries ? 'retry' : 'terminate' }
* Detect whether we're on MASSIVE
    - then set profile
    - unsure if this can be set during runtime, will investigate
* Coding consequences currently not processing alleles in features with compound locations
    - need to better understand limitations of current approach to determine compatibility with these
    - currently not assessing consequences and instead warning user
    - we should be able to handle these if:
        - gene sequence is correctly obtained
        - within SNP position adjusted as needed
            - this will depend on position with in gene
            - likely the most difficult part
        - ensure that the SNP position is not between compound positions
            - currently we only test outer bounds of genes
            - will need to test each bound individual until we're sure it is genic
    - another consideration is whether these annotations have biological meaning
        - do these genes truly have non-contiguous open reading frames?
        - or is this some product of annotation weirdness
* Coding consequences currently only examines each SNP individually
    - report coding result when multiple SNPs are in one codon
    - this probably should replace looking at SNPs individually
    - would require reworking a large amount of `determine_coding_consequences.py`
* Scaling job resources on file size may be more appropriate
    - use the `size()` method, returns file size in bytes
* Test on large dataset
    - haemophilus influenze, currently running
    - typhi, zoe to provide data
* Set output files to be copied or moved rather than symlinked
    - move is preferable, particularly for large runs producing many BAMs
    - could do this through netflow publishDir interface
        - if moving files this will cause problems unless they're not used downstream
    - or through a final stage to run some BASH to move files
        - for file in $(find . -type l); do mv $(readlink ${file}) ${file}; done
* Automate comparison of test output data
* Coding consequences for hets
    - create matrix of hets in the same way we do for homs
    - this should be done conditionally
    - pass the data forward
* Conditionally execute phylogeny process
    - on the basis of having n or more isolates
* Check we have input read sets
    - create different channel to check the first item
* Probably can remove mapped flag check in `get_reads_mapped.awk`
    - was previously passing unfiltered bam
* For `get_snp_sites.awk`, check that the input ref and alt allele can be the same (i.e. not use of '.')
* Revert min. base quality to 20 in the site-specific consensus calling of `create_allele_matrix.py`
    - reduced required quality to use with small test dataset
* Add optional fastqc
    - default behaviour to generate reports
    - option to turn off
* Mixed sample detection
    - SNP:het ratio in replicon statistics
    - use zoe's scripts to plot
        - reads het and q30 vcfs
        - plot %mapped to alt, coloured by hets or homs
* Check user inputs and configuration before executing pipeline


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
* Other mixed sample detection methods?


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
