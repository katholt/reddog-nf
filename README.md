# reddog-nf
A nextflow implementation of the microbial variant calling pipeline [`RedDog`](https://github.com/katholt/RedDog). The
central goal of `RedDog` is to identify high-quality SNPs present in a set of isolates using a mapping-based approach.


## Table of contents
* [Table of contents](#table-of-contents)
* [Quickstart](#quickstart)
* [Full example](#full-example)
* [Usage](#usage)
* [Outputs](#outputs)
* [Requirements](#requirements)
* [Tests](#tests)
* [License](#license)


## Quickstart
### General
For the `RedDog` veterans:
```bash
# Clone the reddog-nf repository
git clone https://github.com/scwatts/reddog-nf.git && cd reddog-nf

# Install dependencies and activate conda environment
conda create -c bioconda -c conda-forge -p $(pwd -P)/conda_env --yes --file config/conda_dependencies.txt
conda activate $(pwd -P)/conda_env

# Set input reads, reference, and output directory in configuration file
vim nextflow.config

# Launch pipeline
./reddog.nf
```


### M3 execution
To run the pipeline on M3, there are additional steps that need to be taken and some important considerations. Software
provision is similar to the general example above but you must first load `miniconda` into your path:
```
# Clone the reddog-nf repository
git clone https://github.com/scwatts/reddog-nf.git && cd reddog-nf

# Load the miniconda3 module and install dependencies
module load miniconda3/4.1.11-python3.5
conda create -c bioconda -c conda-forge -p $(pwd -P)/conda_env --yes --file config/conda_dependencies.txt
```
Note that this is the only version of `miniconda` that will work; other versions are currently broken and will not be fix for
M3 users.

`RedDog` should be run in a `screen` session so that it will persist between ssh sessions:
```
# Create and attach to a screen session, load conda environment
screen -R reddog_example
conda activate $(pwd -P)/conda_env
```

Configuration for input reads, reference, and output directory is done as usual. Once the pipeline has been configured you
can execute the run, however the `-profile` argument must be provided explicitly. This is to prevent inadvertently running
pipeline processes on the head node.
```
# Set input reads, reference, and output directory in configuration file
vim nextflow.config

# Launch pipeline with a specific execution profile
./reddog.nf -profile massive
```


Lastly, the ability to resume pipeline execution is built into nextflow. To resume a run the `-resume` argument is used.
Note that if the output directory already exists, you must confirm that `RedDog` can overwrite its contents by specifying
`--force`.
```
# Reattached to the screen session and resume pipeline execution
screen -R reddog_example
./reddog.nf -profile massive -resume --force
```


## Usage
Running `RedDog` is a two step process that involves configuration and then execution. The `nextflow.config` file controls
pipeline configuration, and execution is done by issuing `./reddog.nf`.

Further detail regarding configuration and execution is given below.


### Configuration
#### Input and output
Running `RedDog` requires a set of short reads in FASTQ format, a reference assembly in GenBank format, and an output
directory. These are set in the **Input and output** section of the `nextflow.config` file, for example:
```bash
// Input and output
reads = 'reads/*.fastq.gz'
reference = 'data/reference.gbk'
output_dir = 'output/'
```

In the above snippet, the input reads are specified with a glob containing the star wildcard character. This expression
selects all files in the `reads/` directory that have a suffix of `.fastq.gz`. If your input files are present in more than
one directory, a list of globs can be provided:
```bash
reads = ['read_set_1/*.fastq.gz', 'read_set_2/*fastq.gz']
```

The nextflow implementation of `RedDog` also handles input containing both paired-end and single-end read sets.


#### Options and parameters
*General*
| Option                            | Description                                                       | Default           |
| ----                              | ----                                                              | ----              |
| `bt2_max_frag_len`                | Maximum fragment length that Bowtie 2 will consider               | 200               |
| `bt2_mode`                        | Bowtie2 preset                                                    | sensitive-local   |
| `var_depth_min`                   | Minimum read depth for vcftuils.pl varFilter to call variants     | 5                 |
| `mapping_cover_min`               | Minimum coverage to pass replicon                                 | 50                |
| `mapping_depth_min`               | Minimum depth to pass replicon                                    | 10                |
| `mapping_mapped_min`              | Minimum percentage mapped to pass replicon, largest replicon only | 50                |
| `outgroup_mod`                    | Modifier for ingroup/outgroup designation                         | 2                 |
| `allele_matrix_support`           | Minimum ratio of read support to call a SNP in allele matrix      | 90                |
| `allele_matrix_cons`              | Minimum % of isolates containing a SNP for it to pass filtering   | 95                |

*HPC*
| Option                            | Description                               | Default   |
| ----                              | ----                                      | ----      |
| `max_retries`                     | Number of times to resubmit job, usually with more resources  | 3         |
| `queue_size`                      | Number of jobs to submit to job scheduler                     | 200       |
| `slurm_account`                   | SLURM account used for job submission                         | js66      |


#### Merge run
A merge run allows you to add additional isolates to an existing `RedDog` output. This is preferable over running an entirely
new `RedDog` analysis as a merge run does not require existing isolate data to be remapped, decreasing runtime considerably.
To execute a merge run you first need to enable merge mode and provide the existing `RedDog` output directory in the **Merge
run** section of the `nextflow.config`, for example:

```bash
// Merge run settings
merge_run = true
existing_run_dir = 'existing_output/'
```

Following this, the new isolate read sets and output directory must be set in the **Input and output** configuration section:
```bash
// Input and output
reads = 'new_reads/*.fastq.gz'
reference = 'data/reference.gbk'
output_dir = 'merged_output/'
```

It is critical that the same reference and filtering/quality settings are used both the previous and new run. There are a
series of preflight checks performed that compare settings and supplied reference between runs. If any inconsistency is
detected, these will be indicated and the pipeline will abort. These inconsistencies can be ignored by setting
`merge_ignore_errors = true`.

At the conclusion of a merge run the output directory (`merged_output/` in the above snippet) will contain all data from the
previous run (from `existing_output/`) and data generated from the new read sets. Once the contents of the new output files
have been checked, the old output directory can be deleted.


### Pipeline execution
#### General
It is recommended that `RedDog` is launched inside a persistent terminal session such as `screen` as the pipeline may run
over several hours or days. In a normal setting without a persistent terminal, closing the terminal or ssh session running
`RedDog` would cause the pipeline to also quit. Executing `RedDog` in a `screen` session allows you to close the terminal or
ssh session and have the pipeline continue in the background. Here is a brief example:
```bash
# Open a new screen session and launch the pipeline
screen -R reddog_run
./reddog.nf

# Detach from the screen session with C-a (ctrl+a)
# List active screen sessions
screen -ls

# Reattach to the screen
screen -R reddog_run
```

Pipeline resume is available through the nextflow framework and enabled by specifying `-resume` on the command line. To avoid
potentially overwriting existing analysis, `RedDog` will refuse to write to an output directory if it already exists unless
you explicitly allow it to do so with `--force`. To resume a run you must provide both arguments:
```bash
./reddog.nf -resume --force
```

The nextflow framework stores intermediate files in a directory called `work/`. This directory can become very large and
should be deleted once the `RedDog` run has completed.


#### M3 specific
When running `RedDog` on M3 you must explicitly set the executor. In the nextflow framework, the executor manages how jobs are
run and the default is to use the `standard` executor, which launches jobs locally. However on M3 jobs should never be run
locally and instead submitted to the SLURM queue. This is done in `RedDog` by using the `massive` executor:
```bash
./reddog.nf -profile massive
```

The pipeline is configured to dynamicly select the optimal partition and QoS for each submitted job. Selection is conditioned
on the walltime of each job:
| Walltime          | QoS       | Partition             |
| ----              | ----      | ----                  |
| 0-30 minutes      | genomics  | comp, genomics, short |
| 31-240 minutes    | genomics  | comp, genomics        |
| 241+ minutes      | normal    | comp                  |

Additionally, resources requested for each task are defined in `config/slurm_job.config`. By default, when a job fails it is
resubmitted two more times but with increased resources. Where a job fails three times due to insufficient resources, you can
manually specify the required resources in `config/slurm_job.config`.

The `genomics` queue limits the number of concurrently queued or running jobs and once this limit is reached, any further job
that is submitted is rejected. To avoid job rejected in this context, the pipeline will only submit 200 tasks to the queue at
any one time. This can be adjusted through the `queue_size` option in `nextflow.config`.


## Outputs
### Overview
For each replicon in the provided reference, the following outputs are generated:
| File                              | Name                                                      |
| ----                              | ----                                                      |
| Mapping statistics                | `reference-name`\_`replicon`\_mapping\_stats.tsv          |
| SNP alleles                       | `reference-name`\_`replicon`\_alleles.csv                 |
| Conserved SNP alleles             | `reference-name`\_`replicon`\_alleles\_cons0.95.csv       |
| Alignment of conserved SNPs       | `reference-name`\_`replicon`\_cons0.95.mfasta             |
| Phylogeny of conserved SNPs       | `reference-name`\_`replicon`\_cons0.95.tree               |
| Consequences of conserved SNPs    | `reference-name`\_`replicon`\_consequences\_cons0.95.tsv  |
| Gene coverage                     | `reference-name`\_gene\_coverage.csv                      |
| Gene depth                        | `reference-name`\_gene\_depth.csv                         |

Additionally several directory outputs are created:
| Directory                         | Location and name         |
| ----                              | ----                      |
| Filtered isolate BAMs             | `output_dir`/bams/        |
| Filtered isolate VCFs             | `output_dir`/vcfs/        |
| Run information                   | `output_dir`/run\_info/   |
| Read quality reports (optional)   | `output_dir`/fastqc/      |


### File descriptions
#### Mapping statistics
| Column name               | Description                                               |
| ----                      | ----                                                      |
| isolate                   | Isolate name                                              |
| replicon\_coverage        | Coverage of replicon for positions >=1 depth              |
| replicon\_average\_depth  | Average depth for positions >=1 depth                     |
| replicon\_reads\_mapped   | Proportion of all mapped reads belonging to this replicon |
| total\_mapped\_reads      | Proportion of reads mapped to any replicon                |
| total\_reads              | Count of total mapped reads                               |
| snps                      | Number of homozygous SNPs                                 |
| snps\_heterozygous        | Number of heterozygous SNPs                               |
| indels                    | Number of indels                                          |
| pass\_fail                | Pass or fail classification based on mapping statistics   |
| phylogeny\_group          | Ingroup or outgroup designation                           |

#### Allele matrix
| Column name       | Description                       |
| ----              | ----                              |
| Pos               | Nucleotide position in reference  |
| Reference         | Reference nucleotide (+ve strand) |
| `isolate columns` | Isolate columns                   |

#### Coding consequences
| Column name                   | Description                               |
| ----                          | ----                                      |
| position                      | Nucleotide position in reference          |
| ref                           | Reference nucleotide (+ve strand)         |
| alt                           | Alternative nucleotide (+ve strand)       |
| strand                        | Encoding strand of gene                   |
| change\_type                  | Mutation type                             |
| gene                          | Gene name                                 |
| ref\_codon                    | Reference codon                           |
| alt\_codon                    | Alternative codon                         |
| ref\_aa                       | Reference amino acid                      |
| alt\_aa                       | Alternative amino acid                    |
| gene\_product                 | Description of encoded protein            |
| gene\_nucleotide\_position    | Nucleotide position of mutation in gene   |
| gene\_codon\_position         | Position of affected codon in gene        |
| codon\_nucleotide\_position   | Nucleotide position of mutation in codon  |
| notes                         | Misc. notes                               |

#### Depth and coverage
| Column name       | Description               |
| ----              | ----                      |
| replicon          | Name of replicon          |
| locus\_tag        | Gene locus tag            |
| `isolate columns` | Isolate columns           |


## Requirements
The following software are required:
* `bcftools`, version ≥1.9
* `biopython`, version ≥1.76
* `bowtie2`, version ≥2.4.1
* `fastqc`, version ≥0.11.9
* `fasttree`, version ≥2.1.10
* `multiqc`, version ≥1.7
* `nextflow`, version ≥20.10.0
* `openjdk`, version ≥8.0.192
* `samtools`, version ≥1.9
* `seqtk`, version ≥1.3


The easiest way to satisfy these dependencies is to use `conda` and the list of dependencies in
`config/conda_dependencies.txt`:
```bash
# First ensure you're in the top level directory of the git repo
conda create -c bioconda -c conda-forge -p $(pwd -P)/conda_env --yes --file config/conda_dependencies.txt
conda activate $(pwd -P)/conda_env
```


## Tests
I've written both end-to-end tests and unit/functional tests. These are designed to detect regressions, test effect of
modifications, and check compatibility of installed software. To run these tests first provision requirements as shown in
[Requirements](#requirements).

```bash
# Unit tests and functional tests, run from the git repo top level directory
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


## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
