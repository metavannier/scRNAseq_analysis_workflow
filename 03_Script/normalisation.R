#-------------------------------------
# Normalisation of the dataset seurat object to annotate
#-------------------------------------

#-------------------------------------
# Library
#-------------------------------------

library(SingleCellExperiment)
library(scran)
library(scuttle)
library(scater)
library(Seurat)
library(Matrix)
library(irlba)
library(data.table)

#-------------------------------------
# Path to files / folders
#-------------------------------------

DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")

SAMPLE_ID = snakemake@params[["sample_id"]]
TEXT_OUTPUT = snakemake@output[["normalisation_output"]]
STEP2 = "02_seurat/"

#----------------------------------------
# Normalisation with Scran
#----------------------------------------

len <- length(SAMPLE_ID)
i = 1

while(i <= len)
{
    so <- readRDS(file.path( OUTPUTDIR, STEP2, SAMPLE_ID[i], paste0(SAMPLE_ID[i],"_filtered_seurat_object.rds")))
    # convert to SingleCellExperiment
    sce_dataset <- as.SingleCellExperiment(so)

    clusts_sce_dataset <- scran::quickCluster(sce_dataset,
                            method = "igraph",
                            use.ranks = FALSE, # suggested by the authors
                            min.size = 100) # require at least 100 cells per cluster

    min.clust = min(table(clusts_sce_dataset))/2
    new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
    sce_dataset = computeSumFactors(sce_dataset, clusters = clusts_sce_dataset, sizes = new_sizes, max.cluster.size = 3000)

    sce_dataset <- logNormCounts(sce_dataset)

    print(rownames(sce_dataset)[1:5])

    print(colData(sce_dataset)[1:5, ])

    so[["RNA"]] <- SetAssayData(so[["RNA"]],
                                  slot = "data", 
                                  new.data = logcounts(sce_dataset))

    so$sizeFactors <- sizeFactors(sce_dataset)

    UpdateSeuratObject(so)

    # Save each seurat object
    saveRDS(so, file = file.path( OUTPUTDIR, STEP2, SAMPLE_ID[i], paste0(SAMPLE_ID[i],"_filteredscran_seurat_object.rds")))

    # Save the normalize matrix as csv file
    normalized_matrix <- GetAssayData(object = so[["RNA"]], slot = "data")
    # normalized_matrix <- as.matrix(normalized_matrix)
    # normalized_matrix <- t(normalized_matrix)
    normalized_matrix <- as.data.frame(normalized_matrix)

    fwrite(x = normalized_matrix , file = file.path( OUTPUTDIR, STEP2, SAMPLE_ID[i], paste0( SAMPLE_ID[i], "_scran_normalized_matrix.csv")), row.names = TRUE)

    # Next sample
    i = i + 1
}

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

output_file<-file(TEXT_OUTPUT)
writeLines(c("Normalisation step finished"), output_file)
close(output_file)