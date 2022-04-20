#!/bin/bash
module load userspace/all
module load python3/3.6.3
#module load singularity/3.5.1

#mkdir snakemake_virtenv
virtualenv-3.6 snakemake_virtenv
source snakemake_virtenv/bin/activate
pip3 install snakemake==5.22.1
install numpy deactivate
#virtualenv-3.6 -p python3 snakemake_virtenv

