#..............................................................
# Preparation of the folowing reference for sims :
# Mouse gastrulation reference :
# doi:10.1038/s41586-019-0933-9
#..............................................................

#-------------------------------------
# Library
#-------------------------------------

library(data.table)
library(SingleCellExperiment)
BiocManager::install("MouseGastrulationData")
library(MouseGastrulationData)
library(scran)
library(scuttle)
library(scater)
library(Seurat)
library(Matrix)
library(irlba)
library(zellkonverter)
library(singleCellTK)

sessionInfo()

#-------------------------------------
# Path to files / folders
#-------------------------------------

DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")

SAMPLE_ID = snakemake@params[["sample_list"]]

MOUSEGASTRULATION_SAMPLES =  snakemake@params[["mousegastrulation_samples"]]

OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]

MOUSEGASTRULATION_MATRIX = snakemake@params[["mousegastrulation_matrix"]]
DATA_MATRIX = snakemake@params[["data_matrix"]]

STEP2 = "02_seurat/"
STEP3 = "03_sims"

#-------------------------------------
# Path to files / folders
#-------------------------------------

##  Loading the samples from the atlas 
MOUSEGASTRULATION_SAMPLES <- type.convert(MOUSEGASTRULATION_SAMPLES, dec=".", as.is = TRUE)
sce <- EmbryoAtlasData(samples = c(MOUSEGASTRULATION_SAMPLES))

# Change ENSEMBL name of the count matrix by corresponding gene name 
sce <- setRowNames(sce, "SYMBOL")

# Remove the cell with Na as labeltype
sce <- sce[,!is.na(colData(sce)$celltype)]

#-------------------------------------
# Normalisation : we follow the process from Pijuan-Sala et al. 2019 
# (doi:10.1038/s41586-019-0933-9)
#-------------------------------------

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

#-------------------------------------
# Keep genes that match both the 
# matrix to annotate and the maytrix of reference
#-------------------------------------

len <- length(SAMPLE_ID)
i = 1
while(i <= len)
{
    # Read genes files from the matrix to annotate
    count_data <- read.delim(file.path( OUTPUTDIR, STEP2, SAMPLE_ID[i], paste0( SAMPLE_ID[i], "_scran_normalized_matrix.csv")), header=TRUE, row.names=1,sep = ",", check.names=FALSE)
    sce_data <- SingleCellExperiment(assays = list(logcounts = count_data))
    rm(count_data)

    # Subset the genes of the matrix to annotate
    keep_genes <- intersect(rownames(sce), rownames(sce_data))
    sims_count_atlas <- sce[match(keep_genes,rownames(sce)), ]
    sims_count_data <- sce_data[na.omit(match(keep_genes,rownames(sce_data))), ]
    rm(sce_data)

    # Transpose rowdata (cells) and column (genes) for sims
    sims_count_atlas_matrix <- as.matrix(logcounts(sims_count_atlas))
    sims_count_atlas_matrix <- t(sims_count_atlas_matrix)
    sims_count_atlas_matrix <- as.data.frame(sims_count_atlas_matrix)

    # Save the atlas normalize matrix as csv file    
    dir.create(file.path(OUTPUTDIR, STEP3, SAMPLE_ID[i]))
    fwrite(x = sims_count_atlas_matrix , file = file.path( OUTPUTDIR, STEP3, SAMPLE_ID[i], paste0( OUTPUT_NAME_REF_MATRIX,".csv")), row.names = TRUE)
    rm(sims_count_atlas_matrix)

    #-------------------------------------
    # Get the label of each cells of the reference atlas
    #-------------------------------------

    coldata <- DataFrame(Cell=colData(sce)$cell,CellType=colData(sce)$celltype)
    coldata <- as.data.frame(coldata)

    fwrite(x = coldata , file = file.path( OUTPUTDIR, STEP3, SAMPLE_ID[i], paste0( OUTPUT_NAME_REF_METADATA,".csv")), row.names = TRUE)
    rm(coldata)

    # For the data to annotate
    sims_count_data <- as.matrix(logcounts(sims_count_data))
    sims_count_data_matrix <- t(sims_count_data)
    sims_count_data_matrix <- as.data.frame(sims_count_data_matrix)
    fwrite(x = sims_count_data_matrix , file = file.path( OUTPUTDIR, STEP3, SAMPLE_ID[i], paste0( OUTPUT_NAME_MATRIX,".csv")), row.names = TRUE)
    rm(sims_count_data)

    # Next sample
    i = i + 1
}

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

output_file<-file(TEXT_OUTPUT)
writeLines(c("Files preparation for the Mouse Gastrulation reference for SIMS finished (CSV format)"), output_file)
close(output_file)
