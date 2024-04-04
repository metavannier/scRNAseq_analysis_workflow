#!/bin/sh

#SBATCH -J allen_SIMS
#SBATCH -p kepler
#SBATCH --ntasks-per-node 24
#SBATCH -A b324
#SBATCH -t 04:00:00
#SBATCH -o ./%N.%x.out
#SBATCH -e ./%N.%x.errs

# This script needs to be started from the run directory

# Load the modules and start the virtual environment
module purge
module load userspace/all
module load python3/3.6.3
module load singularity/3.5.1
pip install snakemake==6.3.0
# module load cuda/11.6

# ================================================
# Run the workflow
# ================================================

# For the worflow in general
snakemake --unlock
snakemake --snakefile Snakefile --use-singularity --use-conda --conda-frontend conda --conda-not-block-search-path-envvars --singularity-args="-B /scratch/$SLURM_JOB_USER/scRNAseq_analysis_workflow/" --cores 24

### For sims
# snakemake --snakefile Snakefile --use-singularity --singularity-args="-B /scratch/$SLURM_JOB_USER/scRNAseq_analysis_workflow/" --cores 18