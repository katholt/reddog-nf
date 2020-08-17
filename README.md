# reddog-nf
A nextflow reimplementation of the microbial variant calling pipeline [RedDog](https://github.com/katholt/RedDog). A
diagrammatic overview of the pipeline can be seen [here](assets/pipeline_overview.pdf).


## Table of contents
* [Table of contents](#table-of-contents)
* [Quickstart](#quickstart)
* [Tests](#tests)
* [Outputs](#outputs)
* [License](#license)


## Quickstart
```bash
# Clone the reddog-nf repository
git clone https://github.com/scwatts/reddog-nf.git && cd reddog-nf

# Install dependencies and activate conda environment
conda create -c bioconda -c conda-forge -p $(pwd -P)/conda_env --yes --file config/conda_dependencies.txt
conda activate $(pwd -P)/conda_env

# Set input reads, reference, and output in configuration file
vim nextflow.config

# Run
./reddog.nf
```


## Tests
First provision the reddog-nf repository and dependencies as shown in [Quickstart](#quickstart).
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


## Outputs
For each replicon in the provided reference, the following outputs are generated:
| File                          | Name                                                  |
| ----                          | ----                                                  |
| Mapping statistics            | `reference_name`\_`replicon`\_mapping\_stats.tsv      |
| All SNP alleles               | `reference_name`\_`replicon`\_alleles.tsv             |
| Conserved SNP alleles         | `reference_name`\_`replicon`\_alleles\_cons.tsv       |
| Conserved SNP alignment       | `reference_name`\_`replicon`\_cons.mfasta             |
| Conserved SNP phylogeny       | `reference_name`\_`replicon`\_cons.tree               |
| Conserved SNP consequences    | `reference_name`\_`replicon`\_consequences\_cons.tsv  |
| Gene depth                    | `reference_name`\_gene\_depth.tsv                     |
| Gene coverage                 | `reference_name`\_gene\_coverage.tsv                  |

Additionally several directory outputs are created:
| Directory                         | Location and name         |
| ----                              | ----                      |
| Filtered isolate BAMs             | `output_dir`/bams/        |
| Filtered isolate VCFs             | `output_dir`/vcfs/        |
| Run information                   | `output_dir`/run\_info/   |
| Read quality reports (optional)   | `output_dir`/fastqc/      |


## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
