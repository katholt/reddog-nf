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


## Full example
* full worked example on M3
* some example data?


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
General:
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

HPC:
| Option                            | Description                               | Default   |
| ----                              | ----                                      | ----      |
| `queue_size`                      | Number of jobs to submit to job scheduler | 200       |
| `slurm_account`                   | SLURM account used for job submission     | js66      |


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

# Deattach from the screen session with C-a (ctrl+a)
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


#### HPC specific
When running `RedDog` on M3 you must explicity set the executor. In the nextflow framework, the executor manages how jobs are
run and the default is to use the `standard` executor, which launches jobs locally. However on M3 jobs should never be run
locally and instead submitted to the SLURM queue. This is done in `RedDog` by using the `massive` executor:
```bash
./reddog.nf --profile massive
```

* dynamic selection of partition and qos
* job submission rate limits
* job resource requests
* SLURM account
* queue size


## Outputs
### Overview
For each replicon in the provided reference, the following outputs are generated:
| File                              | Name                                                      |
| ----                              | ----                                                      |
| Mapping statistics                | `reference_name`\_`replicon`\_mapping\_stats.tsv          |
| SNP alleles                       | `reference_name`\_`replicon`\_alleles.tsv                 |
| Conserved SNP alleles             | `reference_name`\_`replicon`\_alleles\_cons0.95.tsv       |
| Alignment of conserved SNPs       | `reference_name`\_`replicon`\_cons0.95.mfasta             |
| Phylogeny of conserved SNPs       | `reference_name`\_`replicon`\_cons0.95.tree               |
| Consequences of conserved SNPs    | `reference_name`\_`replicon`\_consequences\_cons0.95.tsv  |
| Gene coverage                     | `reference_name`\_gene\_coverage.tsv                      |
| Gene depth                        | `reference_name`\_gene\_depth.tsv                         |

Additionally several directory outputs are created:
| Directory                         | Location and name         |
| ----                              | ----                      |
| Filtered isolate BAMs             | `output_dir`/bams/        |
| Filtered isolate VCFs             | `output_dir`/vcfs/        |
| Run information                   | `output_dir`/run\_info/   |
| Read quality reports (optional)   | `output_dir`/fastqc/      |


### File descriptions
* Described general contents (maybe uses) and columns for each file type

#### Mapping statistics

#### Allele matrix

#### Coding consequences

#### Depth depth and coverage


## Requirements
The following software are required:
* `bcftools`, version ≥1.9
* `biopython`, version ≥1.76
* `bowtie2`, version ≥2.4.1
* `bwa`, version ≥0.7.17
* `ea-utils`, version ≥1.1.2.779
* `fastqc`, version ≥0.11.9
* `fasttree`, version ≥2.1.10
* `multiqc`, version ≥1.7
* `nextflow`, version ≥20.04.1
* `openjdk`, version ≥8.0.192
* `samtools`, version ≥1.9
* `seqtk`, version ≥1.3


This easiest way to satisfy these dependencies is to use `conda` and the list of dependencies in
`config/conda_dependencies.txt`:
```bash
# First ensure you're in the top level directory of the git repo
conda create -c bioconda -c conda-forge -p $(pwd -P)/conda_env --yes --file config/conda_dependencies.txt
conda activate $(pwd -P)/conda_env
```


## Tests
I've written both end-to-end tests and unit/functional tests. These are designed to detect regressions, test effect of
modifications, and check compatibility of install software. To run these tests first provision requirements as shown in
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
