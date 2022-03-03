#!/bin/sh
#SBATCH -J Job_scrna
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A a272
#SBATCH -t 2:00:00
#SBATCH -o ./%N.%x.out
#SBATCH -e ./%N.%x.err

# This script needs to be started from
# the run directory

# Load the modules and start the virtual environment
module load userspace/all
module load python3/3.6.3
module load singularity/3.5.1

# Create a snakemake virtual environment if necessary
if [ ! -d snakemake_virtenv/ ]
then
  bash 04_Workflow/create_snakemake_virtualenv.sh
fi

source snakemake_virtenv/.venv/bin/activate
pip install snakemake
pip install manager
pip install pandas

# move to the working directory
#cd /scratch/$SLURM_JOB_USER/sc-rnaseq/


# ================================================
# Run the workflow
# ================================================

# Run snakemake
snakemake --snakefile Snakefile \
	--reason \
	--use-singularity \
	--singularity-args="-B /scratch/$SLURM_JOB_USER/sc-rnaseq/" \
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

deactivate

