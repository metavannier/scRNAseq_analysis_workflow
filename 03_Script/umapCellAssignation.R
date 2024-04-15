#-------------------------------------
# This script add the cell labels in a seurat object
# Unknown cell labels are add by the use of threshold 
# (difference btw first pred and second pred) 
# UMAP are produced :
# - For each labels
# - Global UMAP with sims annotation without unknown threshold
# - Global UMAP with sims annotation with unknown threshold
#-------------------------------------

library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
RDS_ASSIGNATION = snakemake@params[["rds_assignation"]]
PRED_FILTERED = snakemake@params[["pred_filtered"]]

SAMPLE_ID = snakemake@params[["sample_id"]]

MATRIX = snakemake@params[["matrix"]]

STEP2 = "02_seurat/"
STEP3 = "03_sims/"

# Load data
data <- readRDS(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID,"_filtered_seurat_object.rds")))
pred <- fread(file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_", MATRIX, "_prediction.csv")))

# Assigning predictions and probabilities to data
data$labels <- pred$first_pred
Idents(data) <- data$labels
diff_prob <- pred$first_prob - pred$second_prob
data$diff_prob <- diff_prob

# Threshold to say if a cell is unknown or not
threshold <- snakemake@params[["threshold"]]
cells_threshold <- colnames(data)[data$diff_prob < threshold]

Idents(object = data, cells = cells_threshold) <- "unknown"

# List of unique identifiers
list <- unique(Idents(data))
# list <- c(list,"unknown")

# Function for creating graphics
createUMAPPlot <- function(data, id) {
  cells <- WhichCells(data, idents = id)
  umap <- DimPlot(data, reduction = "umap", group.by = "labels", 
                  cells.highlight = list(cells), cols.highlight = c("brown2"), 
                  cols = "grey") + ggtitle(id)
  return(umap)
}

# Create a list of graphics
plots_list <- lapply(list, createUMAPPlot, data = data)

# Set number of graphics per page
num_graphiques_par_page <- 1

# Divide graphics into pages
pages <- split(plots_list, ceiling(seq_along(plots_list) / num_graphiques_par_page))

#-------------------------------------
# Save each figure on a pdf
#-------------------------------------
pdf(file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_umap_per_labels.pdf")))
for (page in pages) {
  current_page <- plot_grid(plotlist = page, ncol = 1)
  print(current_page)
}
dev.off()

# Create PDF with global UMAP graph (with sims assignation)
title <- paste("UMAP", SAMPLE_ID, sep=' ')
pdf(file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_umap_sims.pdf")), width=12, height=12)
UMAP <- DimPlot(data, group.by = "labels", label.size = 6,repel = TRUE) + theme(legend.position = "bottom") + ggtitle(title)
print(UMAP)
dev.off()

# Create PDF with UMAP graph of cells exceeding threshold (without NA cells)
cells_threshold <- subset(data, subset = diff_prob >= threshold)
title_threshold <- paste("UMAP filtered", SAMPLE_ID, sep=' ')

pdf(file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_umap_sims_threshold.pdf")), width=12, height=12)
UMAP_threshold <- DimPlot(cells_threshold, group.by = "labels", repel = TRUE) + theme(legend.position = "bottom") + ggtitle(title_threshold)
print(UMAP_threshold)
dev.off()

#-------------------------------------
# Update the label in the metadata with
# unknown cells and save the RDS
#-------------------------------------

data$labels <- Idents(data)
saveRDS(data, file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_", RDS_ASSIGNATION)))

#-------------------------------------
# create the predicted matrix with the
# output file for the snakemake rule
#-------------------------------------

tab_pred_filtered <- paste(Cells(cells_threshold),cells_threshold$labels,cells_threshold$diff_prob)
write.table(file=file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_", PRED_FILTERED)), tab_pred_filtered ,sep=",", quote=F, row.names=F, col.names=F)

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

TEXT_OUTPUT <- snakemake@output[["umapAssignation_output"]]

output_file<-file(TEXT_OUTPUT)
writeLines(c("Rules UMAP cell assignation finished"), output_file)
close(output_file)