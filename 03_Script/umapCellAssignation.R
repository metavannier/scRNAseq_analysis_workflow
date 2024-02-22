#..............................................................
# 
#..............................................................

#-------------------------------------
# Library
#-------------------------------------
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]
OUTPUTDIR = file.path((DIRECTORY), "05_Output")

SAMPLE_ID = snakemake@params[["sample_id"]]

MATRIX = snakemake@params[["matrix"]]

STEP2 = "02_seurat/"
STEP3 = "03_sims/"

#-------------------------------------
# Load  the perdiction file and the 
# seurat file
#-------------------------------------
data <- readRDS(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID,"_filtered_seurat_object.rds")))
pred <- fread(file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(MATRIX,"_prediction.csv")))

#-------------------------------------
# Add the first prediction in the 
# seurat object
#-------------------------------------
data$labels <- pred$first_pred
Idents(data) <- data$labels

#-------------------------------------
# Create a list for each label
#-------------------------------------
liste <- unique(pred$first_pred)

#-------------------------------------
# Create a list to stock each UMAP
#-------------------------------------
plots_list <- list()

#-------------------------------------
# To create a UMAP for each label
#-------------------------------------
for (i in liste) {
  cells <- WhichCells(data, idents = c(i))
  umap <- DimPlot(data, reduction = "umap", group.by = "labels", cells.highlight= list(cells), cols.highlight = c("lightblue"), cols= "grey") + ggtitle(i)
  plots_list[[i]] <- umap
}

#-------------------------------------
# Save each figure on a pdf
#-------------------------------------
num_graphiques_par_page <- 1
pages <- split(plots_list, ceiling(seq_along(plots_list) / num_graphiques_par_page))
pdf(file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, "umap_per_labels.pdf"))
for (i in seq_along(pages)) {
    for(plot in pages[[i]]){
        print(plot)
    }
}
dev.off()

#-------------------------------------
# Save the globa figure on a pdf
#-------------------------------------
pdf(file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, "umap.pdf"))
UMAP <- DimPlot(data, group.by = "labels", repel = TRUE) + theme(legend.position = "right")
print(UMAP)
dev.off()

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

TEXT_OUTPUT <- snakemake@output[["umapAssignation_output"]]

output_file<-file(TEXT_OUTPUT)
writeLines(c("Rules UMAP cell assignation finished"), output_file)
close(output_file)