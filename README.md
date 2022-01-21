

# Single-Cell RNA-seq Analysis Workflow for 10x GENOMICS data

## Author

Thomas Vannier (@metavannier), https://centuri-livingsystems.org/t-vannier/

## About

This workflow performs a Snakemake pipeline to process 10x single-cell RNAseq data from fastq files to the differential expression marker-gene analysis.
Correction for technical differences between datasets can be included (i.e. batch effect correction) with the integration method during the sctransform process in Seurat to perform comparative scRNA-seq analysis across experimental conditions.

## Usage

You need to install [Singularity](https://github.com/hpcng/singularity/blob/master/INSTALL.md#install-golang) on your computer. This workflow also work in a slurm environment.

Each snakemake rules call a specific conda environment. In this way you can easily change/add tools for each step if necessary. 

### Step 1: Install workflow

You can use this workflow by downloading and extracting the latest release. If you intend to modify and further extend this workflow or want to work under version control, you can fork this repository.

We would be pleased if you use this workflow and participate in its improvement. If you use it in a paper, don't forget to give credits to the author by citing the URL of this repository and, if available, its [DOI](https://).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files and repositories:
- 00_RawData need the fastq file of each run to analyse.
- 01_Reference need the reference genome in fasta and the corresponding gff for the cellranger mapping step.
- [config.yaml](/config.yaml) indicating the parameters to use.
- Comment the [Snakefile](/Snakefile) on the input line not expected for the pipeline.
- Build the singularity image of cellranger v6.0.0 (file to large for a github repository):

`singularity build cellranger.sif docker://litd/docker-cellranger`

`mv cellranger.sif 02_Container/`

### Step 3: Execute workflow

- You need [Singularity v3.5.3](https://github.com/hpcng/singularity/blob/master/INSTALL.md#install-golang) installed on your computer or cluster.

- Load snakemake from a docker container and run the workflow from the root by using these commands:

`singularity run docker://snakemake/snakemake:v6.3.0`

- Then execute the workflow locally via

`snakemake  --use-conda --use-singularity --cores 10`

### Step 4: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

`snakemake --report global_report.html`
