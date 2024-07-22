#========================================================================
# DEG analysis
#========================================================================

#=====================================
# Library
#=====================================
library(Seurat)
library(ggplot2)
library(dplyr)
# library(ggrepel)
library(data.table)

#=====================================
# Path to file and folder
#=====================================
DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")

TIME_POINT = snakemake@params[["time_point"]]

STEP2 = "02_seurat/"
STEP4 = "04_diffexp/"

dir.create(file.path(OUTPUTDIR, STEP4))

#========================================================================
# Load the data
#========================================================================
WT <- readRDS(file = file.path(OUTPUTDIR, STEP2, paste0(TIME_POINT,"_WT/" ,TIME_POINT, "_WT_seurat_labelTransfert.rds"))) 
KO <- readRDS(file = file.path(OUTPUTDIR, STEP2, paste0(TIME_POINT,"_KO/" ,TIME_POINT, "_KO_seurat_labelTransfert.rds"))) 

# Merge the two condition together
merge_data <- merge(WT, y = KO, add.cell.ids = c("WT", "KO"))

### DE analysis
Idents(merge_data) <- "SCT"

#=====================================
# Per condition
#=====================================
merge_data <- PrepSCTFindMarkers(merge_data)
de_results_per_condition <- FindMarkers(merge_data, ident.1 = "WT", ident.2 = "KO", group.by = "condition", logfc.threshold = -Inf)

de_results_per_condition$diffexpressed <- "Not significant"
de_results_per_condition$diffexpressed[de_results_per_condition$avg_log2FC > 0.6 & de_results_per_condition$p_val_adj < 0.05] <- "Upregulated"
de_results_per_condition$diffexpressed[de_results_per_condition$avg_log2FC < -0.6 & de_results_per_condition$p_val_adj < 0.05] <- "Downregulated"
table(de_results_per_condition$diffexpressed)

de_results_per_condition$delabel <- NA
de_results_per_condition$delabel[de_results_per_condition$diffexpressed != "Not significant"] <- rownames(de_results_per_condition)[de_results_per_condition$diffexpressed != "Not significant"]

# When pvalue = 0 for plotting
de_results_per_condition$plot_p_val_adj <- ifelse(de_results_per_condition$p_val_adj == 0, 1e-300, de_results_per_condition$p_val_adj)

# plot adding up all layers we have seen so far
plot <- ggplot(data=de_results_per_condition, aes(x=avg_log2FC, y=-log10(plot_p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# Save the file in one folder
ggsave(path = file.path(OUTPUTDIR, STEP4, TIME_POINT), filename = paste0("/", TIME_POINT, "_volcanoplot_per_condition.pdf"), plot = plot, device = "pdf")
fwrite(x = de_results_per_condition, file = file.path(OUTPUTDIR, STEP4, TIME_POINT, paste0("/", TIME_POINT, "_diffexp_per_condition.csv")), row.names = TRUE) 

#=====================================
# Per cluster
#=====================================
# Function to have minimum 3 cells per label
filter_cell_types <- function(seurat_obj, min_cells = 3) {
  celltype_counts <- table(seurat_obj$supertype_label)
  valid_celltypes <- names(celltype_counts[celltype_counts >= min_cells])
  seurat_obj_filtered <- subset(seurat_obj, subset = supertype_label %in% valid_celltypes)
  return(seurat_obj_filtered)
}

WT <- filter_cell_types(WT)
KO <- filter_cell_types(KO)

### Keep only the cells type that are in both condition
common_cell_types <- intersect(unique(KO$supertype_label), unique(WT$supertype_label))

KO <- subset(KO, subset = supertype_label %in% common_cell_types)
WT <- subset(WT, subset = supertype_label %in% common_cell_types)

merge_data <- merge(WT, y = KO, add.cell.ids = c("WT", "KO"))

de_list_per_cluster <- list()

for (cell_type in common_cell_types) {
  print(cell_type)
  cell_type_subset <- subset(merge_data, subset = supertype_label == cell_type)
  print(table(cell_type_subset$orig.ident))
  cell_type_subset <- PrepSCTFindMarkers(cell_type_subset)
  de_results <- FindMarkers(cell_type_subset, ident.1 = "WT", ident.2 = "KO", group.by = "condition", logfc.threshold = -Inf)
  de_list_per_cluster[[cell_type]] <- de_results
}

# Assuming de_list_per_cluster is your list containing data frames for each cell type
for (cell_type in names(de_list_per_cluster)) {
  de_results_per_cluster <- de_list_per_cluster[[cell_type]]
  
  # Classify differentially expressed genes
  de_results_per_cluster$diffexpressed <- "Not significant"
  de_results_per_cluster$diffexpressed[de_results_per_cluster$avg_log2FC > 0.6 & de_results_per_cluster$p_val_adj < 0.05] <- "Upregulated"
  de_results_per_cluster$diffexpressed[de_results_per_cluster$avg_log2FC < -0.6 & de_results_per_cluster$p_val_adj < 0.05] <-  "Downregulated"
  
  # Display the table of diffexpressed values
  print(table(de_results_per_cluster$diffexpressed))
  
  # Label differentially expressed genes
  de_results_per_cluster$delabel <- NA
  de_results_per_cluster$delabel[de_results_per_cluster$diffexpressed != "Not significant"] <- rownames(de_results_per_cluster)[de_results_per_cluster$diffexpressed != "Not significant"]
  
  de_results_per_cluster$plot_p_val_adj <- ifelse(de_results_per_cluster$p_val_adj == 0, 1e-300, de_results_per_cluster$p_val_adj)
  
  # Store the modified data frame back in the list
  de_list_per_cluster[[cell_type]] <- de_results_per_cluster
  
  plot <- ggplot(data=de_results_per_cluster, aes(x=avg_log2FC, y=-log10(plot_p_val_adj), col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    # geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") +
    labs(title=paste("Volcano Plot for", cell_type))
  
  # Save the plot
  ggsave(path = file.path(OUTPUTDIR, STEP4, TIME_POINT), filename = paste0("/", TIME_POINT, "_volcanoplot_", cell_type, ".pdf"), plot = plot, device = "pdf")
  fwrite(x = de_results_per_cluster, file = file.path(OUTPUTDIR, STEP4, TIME_POINT, paste0("/", TIME_POINT, "_diffexp_", cell_type, ".csv")), row.names = TRUE) 
}

TEXT_OUTPUT <- snakemake@output[["seurat_DEG_output"]]
output_file<-file(TEXT_OUTPUT)
writeLines(c("Rules seurat DE finished"), output_file)
close(output_file)