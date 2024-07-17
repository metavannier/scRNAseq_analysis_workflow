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
# Loading library
#-------------------------------------

# suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(URD))

# SO <- readRDS(file.path( OUTPUTDIR, STEP, "merged_seurat_object.rds"))
# seuratToURD(SO)