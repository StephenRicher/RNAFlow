# RNA-Flow

## A comprehensive Snakemake workflow for automated processing and quality control of RNA-seq data.

RNA-flow aims to provide an accessible and user-friendly experience to analyse RNA-seq using a wide range of published tools.
The pipeline utilises the workflow management system Snakemake and automatically handles installation of all required software with no user input. RNA-flow can also be easily scaled to work in cluster environments. Current tools utilised by the workflow include:

 * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - A quality control tool for high throughput sequence data.
 * [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) - A tool to screen for species composition in FASTQ sequences.
 * [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) - Remove adapter sequences, primers, poly-A tails and others from high-throughput sequencing data.
 * [Kallisto](https://pachterlab.github.io/kallisto/about) - Quantify abundances of transcripts from RNA-Seq data using pseudoalignment.
 * [HISAT2](https://daehwankimlab.github.io/hisat2/manual/) - A fast and sensitive spliced aligner.
 * [Samtools](http://www.htslib.org/doc/samtools.html) - Interact with high-throughput sequencing data in SAM, BAM and CRAM formats.
 * [deepTools](https://deeptools.readthedocs.io/en/develop/) - Tools for analysis and quality control of high-throughput sequencing data.
 * [RSeQC](http://rseqc.sourceforge.net/) - Comprehensively evaluate and quality check high-throughput RNA sequencing data.
 * [featureCounts](http://subread.sourceforge.net/) - Count reads to genomic features such as genes, exons, promoters and genomic bins.
 * [MultiQC](https://multiqc.info/) - Aggregate results from bioinformatics analyses across many samples into a single report.

## Table of contents

  * [Installation](#installation)
  * [Configuration](#configuration)
  * [Usage](#usage)
     * [Cluster execution](#cluster-execution)
  * [Output](#output)
     * [MultiQC report](#multiqc-report)
     * [Kallisto](#kallisto)
     * [featureCounts](#featureCounts)
  * [References](#references)

## Installation

RNA-Flow works with Python >=3.6 and requires [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

The RNA-Flow repository can be downloaded from GitHub as follows:

```bash
$ git clone https://github.com/StephenRicher/RNA-Flow.git
```

## Configuration
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
Consult the official Snakemake documentation [here](https://snakemake.readthedocs.io/en/latest/executing/cli.html#profiles) to learn more about running RNA-Flow on your particular cluster environment.

## Output

### MultiQC report
RNA-Flow utilises MultiQC to aggregate the QC and metric report across all samples and all compatible tools used in the pipeline. An example MultiQC report produced by RNA-Flow is shown [here](./README_files/multiqc_report.html).

### Kallisto
Kallisto quantifies transcript abundance using pseudoalignment.
The output is intended for downstream processing and differential expression analysis using [Sleuth](https://pachterlab.github.io/sleuth/about).

### featureCounts
featureCounts quantifies reads mapping to genomic features (e.g. exons) following alignment to a reference genome. The output is commonly used for differential expression analysis using tools such as [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [DESeq2](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

## References
* Köster, J. and Rahmann, S., 2012. Snakemake-a scalable bioinformatics workflow engine. Bioinformatics [Online]. Available from: https://doi.org/10.1093/bioinformatics/bts480.
* Hemphill, W., Rivera, O. and Talbert, M., 2018. RNA-sequencing of Drosophila melanogaster head tissue on high-sugar and high-fat diets. G3: Genes, Genomes, Genetics [Online]. Available from: [https://doi.org/10.1534/g3.117.300397](https://doi.org/10.1534/g3.117.300397).
* Simon Andrews, 2020. Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data. Soil.
* Wingett, S.W. and Andrews, S., 2018. FastQ Screen: A tool for multi-genome mapping and quality control. F1000Research [Online]. Available from: https://doi.org/10.12688/f1000research.15931.2.
* Martin, M., 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal [Online]. Available from: https://doi.org/10.14806/ej.17.1.200.
* Bray, N.L., Pimentel, H., Melsted, P. and Pachter, L., 2016. Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology [Online]. Available from: https://doi.org/10.1038/nbt.3519.
* Kim, D., Langmead, B. and Salzberg, S.L., 2015. HISAT: A fast spliced aligner with low memory requirements. Nature Methods [Online]. Available from: https://doi.org/10.1038/nmeth.3317.
* Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G. and Durbin, R., 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics [Online]. Available from: https://doi.org/10.1093/bioinformatics/btp352.
* Ramírez, F., Dündar, F., Diehl, S., Grüning, B.A. and Manke, T., 2014. DeepTools: A flexible platform for exploring deep-sequencing data. Nucleic Acids Research [Online]. Available from: https://doi.org/10.1093/nar/gku365.
* Wang, L., Wang, S. and Li, W., 2012. RSeQC: Quality control of RNA-seq experiments. Bioinformatics [Online]. Available from: https://doi.org/10.1093/bioinformatics/bts356.
* Liao, Y., Smyth, G.K. and Shi, W., 2014. FeatureCounts: An efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics [Online]. Available from: https://doi.org/10.1093/bioinformatics/btt656.
* Ewels, P., Magnusson, M., Lundin, S. and Käller, M., 2016. MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics [Online]. Available from: https://doi.org/10.1093/bioinformatics/btw354.
* Pimentel, H., Bray, N.L., Puente, S., Melsted, P. and Pachter, L., 2017. Differential analysis of RNA-seq incorporating quantification uncertainty. Nature Methods [Online]. Available from: https://doi.org/10.1038/nmeth.4324.
* Robinson, M.D., McCarthy, D.J. and Smyth, G.K., 2009. edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics [Online]. Available from: https://doi.org/10.1093/bioinformatics/btp616.
* Love, M.I., Huber, W. and Anders, S., 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology [Online]. Available from: https://doi.org/10.1186/s13059-014-0550-8.
