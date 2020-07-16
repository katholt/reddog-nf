# reddog-nf

## Quickstart
```bash
# Clone the reddog-nf repository
git clone https://github.com/scwatts/reddog-nf.git && cd reddog-nf

# Install dependencies and activate conda environment
dependencies='''bcftools=1.9
                biopython=1.76
                bowtie2=2.4.1
                bwa=0.7.17
                ea-utils=1.1.2.779
                fastqc=0.11.9
                fasttree=2.1.10
                multiqc=1.7
                nextflow=20.04.1
                openjdk=8.0.192
                samtools=1.9
                seqtk=1.3'''
conda create -p $(pwd -P)/conda_env --yes ${dependencies}
conda activate $(pwd -P)/conda_env

# Configure
vim nextflow.config

# Run
./reddog.nf
```

## Running tests
First clone repo, install dependencies, and activate environment as shown in [Quickstart](#quickstart)
```bash
# Unit tests and functional tests
python3 -m unittest

# End-to-end test
cd tests/simulated_data/
mkdir -p reads/
# Create reference and simulate reads with known mutations
./scripts/generate_reference.py > data/reference.gbk
./scripts/simulate_reads.py --reference_fp data/reference.gbk --spec_fp data/dataset_specification.tsv --output_dir reads/
# Run tests
../../reddog.nf
# Check outputs
./scripts/compare_data.py --run_dir output/ --spec_fp data/dataset_specification.tsv --test_data_dir data/run_data/
```


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
    - reads from bt2 streamed into awk script and placed in either mapped or unmapped bam
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
    - gene counts in isolate set
    - AllStats


## Important note
* Not currently using a strand bias cutoff
    - this is only used in reddog when the runType is set to PE
    - though runType is only ever set to phylogeny or pangenome
    - it may be that runType should be replaced with readType in rr


## TODO
* Set errorStrategy dynamically
    - use `retry` until we reach maximum allowed, then switch to `ignore`
* Currently copying merge data outside of nextflow execution
    - for very large datasets, the script may appear to hang during this operation
    - we could notify user of the operation w progress?
    - or perhaps create a process just to perform the copy w the 'local' executor
* Use nextflow builtin conda functionality
* Table for comparison of differences
* Rewrite gene coverage/depth process to be more efficient
    - determine if more time or memory was required in JH runs
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
