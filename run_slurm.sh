#!/bin/sh
#SBATCH -J Microbeannotator
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -A a272
#SBATCH -t 11:00:00
#SBATCH -o ./%N.%x.out
#SBATCH -e ./%N.%x.err

# This script needs to be started from the run directory

# Load the modules and start the virtual environment
module load userspace/all
module load python3/3.6.3
module load singularity/3.5.1

pip install snakemake==6.3.0
pip install pandas

# ================================================
# Run the workflow
# ================================================

# Run snakemake
snakemake --snakefile Snakefile \
	--reason \
	--use-singularity \
	--use-conda \
	--conda-frontend conda \
	--singularity-args="-B /scratch/$SLURM_JOB_USER/inmed_dechevigny_scrnaseq/" \
	--jobs 1 \
	--latency-wait 20 \
	--max-jobs-per-second 5 \
	--max-status-checks-per-second 5 \
	--cluster-config cluster_config.json \
	--cluster 'sbatch -A {cluster.project} \
		--job-name {cluster.job-name} \
		--partition {cluster.partition} \
		--time {cluster.time} \
		-N {cluster.nodes-number} \
		-n {cluster.cores} \
		--mem-per-cpu {cluster.mem-per-cpu} \
		--output {cluster.output} \
		--error {cluster.error}' \
