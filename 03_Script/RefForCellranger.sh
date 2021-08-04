#!/bin/bash

# Building the reference transcriptome from https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build

## Load parameters
# Input files
link_ref_fasta=$1
link_ref_gff=$2
# Output files
ref_cellranger=$3
# Output path
path_ref=$4
# Reference name
ref_name=$5
# Reference version
ref_version=$6

# Get the files name
split_link_ref_fasta=(${link_ref_fasta//\// })
fasta_name_gz=${split_link_ref_fasta[-1]}
fasta_name=${fasta_name_gz//".gz"/}

split_link_ref_gff=(${link_ref_gff//\// })
gff_name_gz=${split_link_ref_gff[-1]}
gff_name=${gff_name_gz//".gz"/}

# File name gtf filtered output
gff_name_filtered=${gff_name//".gtf"/".filtered.gtf"}

cd $4

Download and unzip
wget ${link_ref_fasta}
gunzip ${fasta_name_gz}

wget ${link_ref_gff}
gunzip ${gff_name_gz}

cellranger mkgtf ${gff_name} ${gff_name_filtered} \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense

cellranger mkref --genome=$5 \
                 --fasta=${fasta_name} \
                 --genes=${gff_name_filtered} \
                 --ref-version=${ref_version}