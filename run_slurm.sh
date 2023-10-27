#!/bin/sh
#SBATCH -J Snakemake_run
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -A b324
#SBATCH -t 2:00:00
#SBATCH -o ./%N.%x.out
#SBATCH -e ./%N.%x.err

# This script needs to be started from the run directory

# Load the modules and start the virtual environment
module purge
module load userspace/all
module load python3/3.6.3
module load singularity/3.5.1
pip install snakemake==6.3.0
pip install pandas

# ================================================
# Run the workflow
# ================================================

snakemake --snakefile Snakefile --use-singularity --use-conda --conda-frontend conda --conda-not-block-search-path-envvars --singularity-args="-B /scratch/$SLURM_JOB_USER/ibdm_lenne_scrnaseq_gastruloids/" --cores 12