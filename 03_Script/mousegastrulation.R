#..............................................................
# Preparation of the folowing reference for sims :
# Mouse gastrulation reference :
# doi:10.1038/s41586-019-0933-9
#..............................................................

#-------------------------------------
# Library
#-------------------------------------

# library(magick)
# # library(BiocManager)
# # library(reshape2)
# # library(data.table)
# library(Matrix)
# library(dplyr)
library(SingleCellExperiment)
BiocManager::install("MouseGastrulationData")
library(MouseGastrulationData)
library(scran)
library(scuttle)
library(scater)
# library(Matrix)
# library(igraph)

sessionInfo()

#-------------------------------------
# Path to files / folders
#-------------------------------------

DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")

SAMPLE_ID = snakemake@params[["sample_id"]]

MOUSEGASTRULATION_SAMPLES =  snakemake@params[["mousegastrulation_samples"]]

OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]

STEP2 = "02_seurat/"
STEP3 = "03_sims"

#-------------------------------------
# Path to files / folders
#-------------------------------------

##  Loading the samples from the atlas 
MOUSEGASTRULATION_SAMPLES <- type.convert(MOUSEGASTRULATION_SAMPLES, dec=".")

sce <- EmbryoAtlasData(samples = c(MOUSEGASTRULATION_SAMPLES))

#-------------------------------------
# Normalisation : we follow the process from Pijuan-Sala et al. 2019 
# (doi:10.1038/s41586-019-0933-9)
#-------------------------------------


## Calculate size factors for normalisation with scran.
# #### DANS LE GIT HUB, VOIR SI A FAIRE SACHANT QUE LES SIZEFACTOR SONT DANS L'OBJET sce
# lib.sizes = Matrix::colSums(counts(sce))
## calcAverage doesn't work
# sce = sce[calcAverage(sce)>0.1,]

## STEP from https://rstudio-pubs-static.s3.amazonaws.com/699579_c6be4bf3220746088bdfd12a61aa15c4.html
# and https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/4_normalisation/normalise.Rmd
# For pre-clustering, we use scran's `quickCluster` function, using the `igraph` method. We specify a maximum cluster size of 3000 cells and a minimum cluster size of 100 cells.
clusts <- as.numeric(quickCluster(sce,
                        method = "igraph",
                        use.ranks = FALSE, # suggested by the authors
                        min.size = 100)) # require at least 100 cells per cluster
# Number of cells in each cluster should be at least twice that of the largest 'sizes'
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

## Use scuttle to normalize function, which calculates log2-transformed normalized expression values.
# This is done by dividing each count by its size factor, adding a pseudo-count and log-transforming.
sce <- logNormCounts(sce)

logcounts(sce)[40:50, 40:50]

# LA mettre objet seurat des data de Pierre François et voir la différence entre normalisation seurat et scutle


# #-------------------------------------
# # Metadata reference file
# #-------------------------------------

# ## Produce the reference metadata matrice
# reference_metadata <- cbind(sce$cell,sce$celltype)

# #-------------------------------------
# # Wtrite the two new matrix +
# # create directory
# #-------------------------------------

# dir.create(file.path(OUTPUTDIR, STEP3, SAMPLE_ID))

# class(sce) <- "numeric"
# write.csv(normcounts_atlas, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_MATRIX,".csv")))

# fwrite(x = reference_metadata, file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_METADATA,".csv")))

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

output_file<-file(TEXT_OUTPUT)
writeLines(c("Files preparation for the Mouse Gastrulation reference for SIMS finished (CSV format)"), output_file)
close(output_file)
