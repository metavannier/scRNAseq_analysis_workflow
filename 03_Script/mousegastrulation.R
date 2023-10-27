#..............................................................
# Preparation of the folowing reference for sims :
# Mouse gastrulation reference :
# doi:10.1038/s41586-019-0933-9
#..............................................................

#-------------------------------------
# Library
#-------------------------------------

library(magick)
library(BiocManager)

BiocManager::install("MouseGastrulationData")
library(MouseGastrulationData)

#-------------------------------------
# Path to files / folders
#-------------------------------------

TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]

#-------------------------------------
# create the output file
#-------------------------------------

output_file<-file(TEXT_OUTPUT)
writeLines(c("Files preparation for the Mouse Gastrulation reference for SIMS finished (CSV format)"), output_file)
close(output_file)
