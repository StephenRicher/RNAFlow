# RNA-Flow

## A comprehensive Snakemake workflow for automated processing and quality control of RNA-seq data.

RNA-flow aims to provide an accessible and user-friendly experience to analyse RNA-seq using a wide range of published tools.
The pipeline utilises the workflow management system Snakemake and automatically handles installation of all required software with no user input. RNA-flow can also be easily scaled to work in cluster environments. Current tools utilised by the workflow include:

 * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - A quality control tool for high throughput sequence data.
 * [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) - A tool to screen for species composition in FASTQ sequences.
 * [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) - A tool to remove adapter sequences, primers, poly-A tails and others from high-throughput sequencing reads.
 * [MultiQC](https://multiqc.info/) - Aggregate results from bioinformatics analyses across many samples into a single report.

## Table of contents

  * [Installation](#installation)
  * [Configuration](#configuration)
  * [Usage](#usage)
  * [Output](#output)
     * [MultiQC report](#multiqc-report)
  * [References](#references)

## Installation

RNA-Flow works with Python >=3.6 and requires [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

The RNA-Flow repository can be downloaded from GitHub as follows:

```bash
$ git clone https://github.com/StephenRicher/RNA-Flow
```

## Configuring RNA-Flow
The RNA-Flow pipeline is fully controlled through a single configuration file that describes parameter settings and paths to relevant files in the system.
RNA-Flow is bundled with a downsampled RNA-seq dataset from *Drosophila melanogaster* (chr3L) ([Hemphill et al., 2018](https://doi.org/10.1534/g3.117.300397)) to test and serve as a template for configuring other datasets.
The configuration file for this example dataset is shown below and can be found at `example/config/config.yaml`.

**Note:** If relative file paths are provided in the configuration file then these are **relative to the working directory**.
The working directory itself (defined by workdir) is relative to the directory Snakemake is executed in.
If not set, the working directory defaults to the directory containing the Snakefile.
Relative paths can be confusing, they are used here to ensure the example dataset works for all users.
If in doubt simply provide absolute paths.

```bash
# Specify output directory (path relative to Snakefile).
# All other paths in config file are relative to this working directory.
workdir: example/analysis/

# Path to table describing data paths.
data:  ../config/samples.tsv

# Set True if paired sequences provided.
paired: False

# For paired-end stranded set to either - RF or FR
# For single-end stranded set to either - R or F
# For unstranded set to 'unstranded'
strand: R

# Set genome configuration.
genome:
    build: BDGP6
    transcriptome:
        - ../genome/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz
        - ../genome/Drosophila_melanogaster.BDGP6.32.ncrna.fa.gz
    index: ../genome/index/BDGP6_tran # Prefix for hisat2 index files
    gff3: ../genome/Drosophila_melanogaster.BDGP6.32.104.gff3.gz

# Cutadapt parameters - see https://cutadapt.readthedocs.io/en/stable/guide.html
cutadapt:
    forwardAdapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    reverseAdapter: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
    overlap:        3
    errorRate:      0.1
    minimumLength:  1
    qualityCutoff:  0,0
    GCcontent:      43

misc:
    # Downsample to N reads before alignment. Useful if aligning for QC only.
    # Set to 0 for no downsampling.
    sample:                    0
    # Reduce for example dataset. Not usually required.
    inferExperimentSampleSize: 75000

# Custom configuration file for multiQC.
multiQCconfig: ../config/multiqc_config.yaml

# Configuration file of paths to genome indexes for FastQ Screen.
# See template in example/config/fastqScreen.config.
fastqScreen:

# Set directory for temporary file. Currently unused.
tmpdir: example/analysis/
```

## Usage

Once Snakemake is installed the example dataset can be processed using the following command.
This command should be run from the RNA-Flow base directory, containing the Snakefile.

```bash
$ snakemake --use-conda --cores 4 --configfile example/config/config.yaml
```

This command will first install all of the relevant Conda environments within the defined working directory (`example/analysis/`).
This may take some time.
The pipeline should then run to completion producing the same figures as shown in the example output below.
Alternatively, you may also want to install the Conda environments in a custom directory.
This is useful when you are running multiple independent analyses and do not wish to repeatedly install the same Conda environments.

```bash
$ snakemake --use-conda --conda-prefix /path/envs/ --cores 4 --configfile example/config/config.yaml
```

### Cluster Execution
All Snakemake based pipelines, including RNA-Flow, are compatible with cluster environments.
Consult the official Snakemake documentation [here](https://snakemake.readthedocs.io/en/v5.25.0/executing/cli.html#profiles) to learn more about running RNA-Flow on your particular cluster environment.

## Output

### MultiQC report

RNA-Flow utilises MultiQC to aggregate the QC and metric report across all samples and all compatible tools used in the pipeline. An example MultiQC report produced by RNA-Flow is shown [here](./README_files/multiqc_report.html).  


###### References
Hemphill, W., Rivera, O. and Talbert, M., 2018. RNA-sequencing of Drosophila melanogaster head tissue on high-sugar and high-fat diets. G3: Genes, Genomes, Genetics [Online]. Available from: [https://doi.org/10.1534/g3.117.300397](https://doi.org/10.1534/g3.117.300397).
