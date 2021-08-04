

# Single-Cell RNA-seq Analysis Workflow for 10x GENOMICS data

## Author

Thomas Vannier (@metavannier), https://centuri-livingsystems.org/t-vannier/

### TO BE CHANGE

## About

This workflow performs a Snakemake pipeline for cellranger to process 10x single-cell RNAseq data from fastq files to the production of a report with result of the marker-gene analysis. The analysis are mainly been achieved using plugins from [QIIME 2](https://qiime2.org/), [Bolyen et al., 2019](https://doi.org/10.1038/s41587-019-0209-9).

You need to install [Singularity](https://github.com/hpcng/singularity/blob/master/INSTALL.md#install-golang) on your computer. This workflow also work in a slurm environment.

Each snakemake rules call a specific conda environment. In this way you can easily change/add tools for each step if necessary. 

## Usage

### Step 1: Install workflow

You can use this workflow by downloading and extracting the latest release. If you intend to modify and further extend this workflow or want to work under version control, you can fork this repository.

We would be pleased if you use this workflow and participate in its improvement. If you use it in a paper, don't forget to give credits to the author by citing the URL of this repository and, if available, its [DOI](https://doi.org/10.5281/zenodo.4772710).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files and repositories:
- 00_RawData need the fastq file of each run to analyse
- 01_Reference the [Marker gene reference databases](https://docs.qiime2.org/2021.4/data-resources/) and the [SEPP reference databases](https://github.com/smirarab/sepp-refs/) if you build the phylogenetic tree with SEPP.
- [manifest.tsv](/manifest.tsv), [sample-metadata.tsv](/sample-metadata.tsv) and [conditions.tsv](/conditions.tsv) to indicate the samples, run, metadata and conditions for the analyse.
- [config.yaml](/config.yaml) indicating the parameters to use.
- Comment the [Snakefile](/Snakefile) on the input line not expected for the pipeline. This depend of the denoising, phylogenetic methods used and if longitudinal analysis are possible (to be improved).

### Step 3: Execute workflow

- You need [Singularity v3.5.3](https://github.com/hpcng/singularity/blob/master/INSTALL.md#install-golang) installed on your computer or cluster.

- Load snakemake from a docker container and run the workflow from the root by using these commands:

`singularity run docker://snakemake/snakemake:v6.3.0`

- Then execute the workflow locally via

`snakemake  --use-conda --use-singularity --cores 10`

### Step 4: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

`snakemake --report report.html`
