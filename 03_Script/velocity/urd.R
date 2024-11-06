#..............................................................
# URD : URD is an R package for reconstructing transcriptional 
# trajectories underlying specification or differentiation 
# processes in the form of a branching tree, using single 
# cell RNA-sequencing data. URD uses a diffusion map projection 
# of the data and works by simulating biased random walks 
# through that projection. 
#..............................................................

#-------------------------------------
# Path to files / folders
#-------------------------------------

DIRECTORY = getwd()
OUTPUTDIR = "05_Output"
STEP = "04_velocity/"
TEXT_OUTPUT = snakemake@output[["urd_output"]]

#-------------------------------------
# Installing packages 
# when not possible with conda env
#-------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org/"))
install.packages("BiocManager")
BiocManager::install("GenomeInfoDbData")
devtools::install_github("farrellja/URD")
# remotes::install_version("Seurat", "3.2.3")

#-------------------------------------
# Loading library
#-------------------------------------

suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(URD))
# library(Seurat)

knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()

SO <- readRDS(file.path( OUTPUTDIR, STEP, "merged_seurat_object.rds"))
seuratToURD(SO)

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

TEXT_OUTPUT <- snakemake@output[["urd_output"]]

output_file<-file(TEXT_OUTPUT)
writeLines(c("Rules URD finished"), output_file)
close(output_file)