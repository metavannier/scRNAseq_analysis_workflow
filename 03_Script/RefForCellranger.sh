#!/bin/bash

# Building the reference transcriptome from https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build

#----------------------------------------
# Load parameters 
#----------------------------------------

# Input files
link_ref_fasta=$1
link_ref_gtf=$2

# Output files
ref_cellranger=$3

# Output path
path_ref=$4

# Reference name
ref_name=$5

# Reference version
ref_version=$6

#----------------------------------------
# Get the files name 
#----------------------------------------

# Fasta
split_link_ref_fasta=(${link_ref_fasta//\// })
fasta_name_gz=${split_link_ref_fasta[-1]}
fasta_name=${fasta_name_gz//".gz"/}

# Gtf
split_link_ref_gtf=(${link_ref_gtf//\// })
gtf_name_gz=${split_link_ref_gtf[-1]}
gtf_name=${gtf_name_gz//".gz"/}

#----------------------------------------
# File name gtf filtered output 
#----------------------------------------

gtf_name_filtered=${gtf_name//".gtf"/".filtered.gtf"}

#----------------------------------------
# Download ans unzip 
#----------------------------------------

cd $4

# Fasta
wget ${link_ref_fasta}
gunzip ${fasta_name_gz}

# Gtf
wget ${link_ref_gtf}
gunzip ${gtf_name_gz}

#----------------------------------------
# Cell Ranger filter with mkgtf
#----------------------------------------

cellranger mkgtf ${gtf_name} ${gtf_name_filtered} \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense

# rm -r ${ref_name}

#----------------------------------------
# Cell Ranger custom reference with mkref
#----------------------------------------

cellranger mkref --genome=$5 \
                 --fasta=${fasta_name} \
                 --genes=${gtf_name_filtered} \
                 --ref-version=${ref_version}