#========================================================================
# Label transfert seurat
#========================================================================

#=====================================
# Library
#=====================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)

#=====================================
# Path to file and folder
#=====================================
DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
REF = file.path((DIRECTORY), "01_Reference")

SAMPLE_ID = snakemake@params[["sample_id"]]

STEP2 = "02_seurat/"

#========================================================================
# Load the data
#========================================================================

### Reference dataset Allen 
ref_data <- readRDS(file = file.path(REF, "/remi/refA21_all.rds"))
premier_modele <- ref_data@assays[["SCT"]]@SCTModel.list[[1]]
ref_data@assays[["SCT"]]@SCTModel.list <- list(premier_modele)

# ### Reference dataset arlotta
# ref_data <- readRDS(file = file.path(REF, "/modify_arlotta_seurat_object.rds"))
# ref_data <- PercentageFeatureSet(ref_data, pattern = "^mt-", col.name = "percent.mt")
# ref_data <- SCTransform(ref_data, vars.to.regress = "percent.mt", assay = "RNA", new.assay.name = "SCT", verbose = FALSE)

### Querry dataset
querry_data <- readRDS(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "_filtered_seurat_object.rds")))### To change here
querry_data <- SCTransform(querry_data, vars.to.regress = "percent.mt", assay = "RNA", new.assay.name = "SCT", verbose = FALSE)

# Change the default assay to SCT
DefaultAssay(ref_data) <- "SCT"
DefaultAssay(querry_data) <- "SCT"

# Keep genes that are in common between querry and ref datasets
ref_genes <- as.data.frame(rownames(ref_data))
querry_genes <- as.data.frame(rownames(querry_data))
all_genes <- intersect(ref_genes$`rownames(ref_data)`, querry_genes$`rownames(querry_data)`)

#========================================================================
# Create functions
#========================================================================

# For findTransfertAnchors function
calculate_k_score <- function(num_cells, proportion = 0.1, min_value = 1) {
  k_score <- max(min_value, floor(num_cells * proportion))
  return(k_score)
}

# For transferData function
calculate_k_weight <- function(sample_size, proportion = 0.2, min_value = 1) {
  k_weight <- round(sample_size * proportion)
  return(max(min_value, k_weight))
}

#========================================================================
# Step 1 : annotate with broad class
#========================================================================

# Create a new column in metadata
querry_data$class_label <- NA

K_SCORE <- 30 # Seurat option

if(ncol(querry_data) < K_SCORE){

  # Calculate a new k_score
  k_score_reference <- calculate_k_score(num_cells = ncol(ref_data))
  k_score_query <- calculate_k_score(num_cells = ncol(querry_data))
  k_score <- min(k_score_reference, k_score_query)

  anchors <- FindTransferAnchors(reference = ref_data,
                               query = querry_data,
                               k.score = k_score,
                               normalization.method = "SCT",
                               features = all_genes)
} else{
  anchors <- FindTransferAnchors(reference = ref_data,
                               query = querry_data,
                               normalization.method = "SCT",
                               features = all_genes)
}

# Number of anchors
num_anchors <- nrow(anchors@anchors)

if(ncol(querry_data) < num_anchors){

  # Calculate the k weight
  k_weight <- calculate_k_weight(sample_size = num_anchors)

  # Prediction with transfertData
  prediction <- TransferData(anchorset = anchors,
                           refdata = ref_data$class_label,
                           k.weight = k_weight
                           )
} else{
  prediction <- TransferData(anchorset = anchors,
                           refdata = ref_data$class_label)
}

querry_data <- AddMetaData(querry_data, metadata = prediction)
querry_data$class_label[Cells(querry_data)] <- querry_data$predicted.id

# Plot UMAP

pdf(file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "_seurat_labelTransfert.pdf")))
DimPlot(querry_data, group.by = "class_label", label.size = 6) + theme(legend.position = "bottom",legend.text=element_text(size=4))

# Plot UMAP per labels
labels <- unique(querry_data@meta.data$`class_label`)

for (label in labels) {
  plot <- DimPlot(querry_data, group.by = "class_label", cells.highlight = WhichCells(querry_data, expression = class_label == label)) + ggtitle(label)
  print(plot)
}

### UMAP plot with the prediction
FeaturePlot(querry_data, features = "prediction.score.max", reduction = "umap")



#========================================================================
# Step 2 : apply a threshold
#========================================================================

ggplot(querry_data@meta.data, aes(x = prediction.score.max)) + 
  geom_histogram(binwidth = 0.05) + 
  theme_minimal() + 
  xlab("Prediction Score Max") + 
  ylab("Frequency") + 
  ggtitle("Distribution of Prediction Scores (Broad Class)")

dev.off()

# Calculate mean and standard deviation
mean_score <- mean(querry_data@meta.data$prediction.score.max)
print(mean_score)

sd_score <- sd(querry_data@meta.data$prediction.score.max)
print(sd_score)

# Define threshold
threshold <- mean_score - sd_score
print(threshold)

# Filter cells
querry_data <- subset(querry_data, subset = prediction.score.max >= threshold)
gc()

pdf(file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "_seurat_labelTransfert.pdf")))

DimPlot(querry_data, group.by = "class_label", label.size = 6) + theme(legend.position = "bottom",legend.text=element_text(size=4))

labels <- unique(querry_data@meta.data$`class_label`)

for (label in labels) {
  plot <- DimPlot(querry_data, group.by = "class_label", cells.highlight = WhichCells(querry_data, expression = class_label == label)) + ggtitle(label)
  print(plot)
}

dev.off()

#========================================================================
# Step 3 : annotate with supertype
#========================================================================

querry_data$supertype_label <- NA

for (class in unique(querry_data$class_label)) {
  print("The supertype class :")
  print(class)
  
  # Set idents=
  Idents(querry_data) <- "class_label"
  Idents(ref_data) <- "class_label"
  
  # Per class
  subset_querry_data <- subset(querry_data, idents = class)
  subset_ref_data <- subset(ref_data, idents = class)
  
  # Check subset sizes
  print(paste("Number of cells in query subset:", ncol(subset_querry_data)))
  print(paste("Number of cells in reference subset:", ncol(subset_ref_data)))

  # Keep the same genes
  ref_genes <- as.data.frame(rownames(subset_ref_data))
  querry_genes <- as.data.frame(rownames(subset_querry_data))
  all_genes <- intersect(ref_genes$`rownames(subset_ref_data)`, querry_genes$`rownames(subset_querry_data)`)
  
  if(ncol(subset_querry_data) < K_SCORE){

  # Calculate a new k_score
  k_score_reference <- calculate_k_score(num_cells = ncol(subset_ref_data))
  k_score_query <- calculate_k_score(num_cells = ncol(subset_querry_data))
  k_score <- min(k_score_reference, k_score_query)

  anchors <- FindTransferAnchors(reference = subset_ref_data,
                               query = subset_querry_data,
                               k.score = k_score,
                               normalization.method = "SCT",
                               features = all_genes)
} else{
  anchors <- FindTransferAnchors(reference = subset_ref_data,
                               query = subset_querry_data,
                               normalization.method = "SCT",
                               features = all_genes)
}

# Number of anchors
num_anchors <- nrow(anchors@anchors)

if(ncol(subset_querry_data) < num_anchors){

  # Calculate the k weight
  k_weight <- calculate_k_weight(sample_size = num_anchors)

  # Prediction with transfertData
  prediction <- TransferData(anchorset = anchors,
                           refdata = subset_ref_data$supertype_label,
                           k.weight = k_weight
                           )
} else{
  prediction <- TransferData(anchorset = anchors,
                           refdata = subset_ref_data$supertype_label)
}
  
  # Add to the metadata
  subset_querry_data <- AddMetaData(subset_querry_data, metadata = prediction)
  querry_data$supertype_label[Cells(subset_querry_data)] <- subset_querry_data$predicted.id
}
pdf(file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "_seurat_labelTransfert.pdf")))


### Plot UMAP
Idents(querry_data) <- querry_data$supertype_label
DimPlot(querry_data, group.by = "supertype_label", label.size = 6) + theme(legend.position = "bottom",legend.text=element_text(size=4))

### Plot UMAP per labels
labels <- unique(querry_data@meta.data$`supertype_label`)

for (label in labels) {
  plot <- DimPlot(querry_data, group.by = "supertype_label", cells.highlight = WhichCells(querry_data, expression = supertype_label == label)) + ggtitle(label)
  print(plot)
}

### UMAP plot with the prediction
FeaturePlot(querry_data, features = "prediction.score.max", reduction = "umap")

dev.off()

saveRDS(querry_data, file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID,paste0(SAMPLE_ID, "_seurat_labelTransfert.rds")))


TEXT_OUTPUT <- snakemake@output[["seurat_labelTransfert_output"]]

output_file<-file(TEXT_OUTPUT)
writeLines(c("Rules seurat Label Transfert finished"), output_file)
close(output_file)



################################################################################################ LABEL TRANSFERT FROM REMI'S ANNOTATION : WT TO KO

# #========================================================================
# # Load the data
# #========================================================================

# ### Ouvrir la matrice de référence pour le NN
# remi_nn <- readRDS("/scratch/lchabot/inmed_deChevigny_scrnaseq/scRNAseq_analysis_workflow/01_Reference/remi/GEcortex_subRef_cluster_class.list.rds")
# remi_nn <- remi_nn[["P30.scRNA-seq.10X"]]
# gc()

# ### Matrice à annoter 
# # querry_data <- readRDS("/scratch/lchabot/inmed_deChevigny_scrnaseq/scRNAseq_analysis_workflow/05_Output/02_seurat/P5_KO/P5_KO_filtered_seurat_object.rds")
# querry_data <- readRDS("/scratch/lchabot/inmed_deChevigny_scrnaseq/scRNAseq_analysis_workflow/05_Output/02_seurat/P30_KO_integrated/P30_KO_integrated_filtered_seurat_object.rds")

# querry_data <- SCTransform(querry_data, vars.to.regress = "percent.mt", assay = "RNA", new.assay.name = "SCT", verbose = FALSE)
# gc()

# ### Set the default assay to SCT for booth seurat object
# DefaultAssay(remi_nn) <- "SCT"
# DefaultAssay(querry_data) <- "SCT"

# # Keep genes that are in common between querry and ref datasets
# ref_genes <- as.data.frame(rownames(remi_nn))
# querry_genes <- as.data.frame(rownames(querry_data))
# all_genes <- intersect(ref_genes$`rownames(remi_nn)`, querry_genes$`rownames(querry_data)`)

# #========================================================================
# # Create functions
# #========================================================================

# # For findTransfertAnchors function
# calculate_k_score <- function(num_cells, proportion = 0.1, min_value = 1) {
#   k_score <- max(min_value, floor(num_cells * proportion))
#   return(k_score)
# }

# # For transferData function
# calculate_k_weight <- function(sample_size, proportion = 0.2, min_value = 1) {
#   k_weight <- round(sample_size * proportion)
#   return(max(min_value, k_weight))
# }

# #========================================================================
# # Step 1 : annotate broad class 
# #========================================================================

# # Create a new column in metadata
# querry_data$class_label <- NA

# K_SCORE <- 30 # Seurat option

# if(ncol(querry_data) < K_SCORE){

#   # Calculate a new k_score
#   k_score_reference <- calculate_k_score(num_cells = ncol(remi_nn))
#   k_score_query <- calculate_k_score(num_cells = ncol(querry_data))
#   k_score <- min(k_score_reference, k_score_query)

#   anchors <- FindTransferAnchors(reference = remi_nn,
#                                query = querry_data,
#                                k.score = k_score,
#                                normalization.method = "SCT",
#                                features = all_genes)
# } else{
#   anchors <- FindTransferAnchors(reference = remi_nn,
#                                query = querry_data,
#                                normalization.method = "SCT",
#                                features = all_genes)
# }

# # Number of anchors
# num_anchors <- nrow(anchors@anchors)

# if(ncol(querry_data) < num_anchors){

#   # Calculate the k weight
#   k_weight <- calculate_k_weight(sample_size = num_anchors)

#   # Prediction with transfertData
#   prediction <- TransferData(anchorset = anchors,
#                            refdata = remi_nn$class_label,
#                            k.weight = k_weight
#                            )
# } else{
#   prediction <- TransferData(anchorset = anchors,
#                            refdata = remi_nn$class_label)
# }

# querry_data <- AddMetaData(querry_data, metadata = prediction)
# querry_data$class_label[Cells(querry_data)] <- querry_data$predicted.id

# # Plot UMAP
# DimPlot(querry_data, group.by = "class_label", label.size = 6) + theme(legend.position = "bottom",legend.text=element_text(size=4))

# # Plot UMAP per labels
# labels <- unique(querry_data@meta.data$`class_label`)

# for (label in labels) {
#   plot <- DimPlot(querry_data, group.by = "class_label", cells.highlight = WhichCells(querry_data, expression = class_label == label)) + ggtitle(label)
#   print(plot)
# }

# ### UMAP plot with the prediction
# FeaturePlot(querry_data, features = "prediction.score.max", reduction = "umap")


# #=============================================
# # Apply a threshold
# #=============================================

# ### Remove cells that are not in the rest of the studies
# querry_data <- subset(querry_data, subset = class_label != "Non-Neuronal")
# querry_data <- subset(querry_data, subset = class_label != "Immature/Migrating")
# querry_data <- subset(querry_data, subset = class_label != "Undetermined")
# querry_data <- subset(querry_data, subset = class_label != "CR")

# DimPlot(querry_data, group.by = "class_label")

# ggplot(querry_data@meta.data, aes(x = prediction.score.max)) + 
#   geom_histogram(binwidth = 0.05) + 
#   theme_minimal() + 
#   xlab("Prediction Score Max") + 
#   ylab("Frequency") + 
#   ggtitle("Distribution of Prediction Scores (Broad Class)")

# # Calculate mean and standard deviation
# mean_score <- mean(querry_data@meta.data$prediction.score.max)
# print(mean_score)

# sd_score <- sd(querry_data@meta.data$prediction.score.max)
# print(sd_score)

# # Define threshold
# threshold <- mean_score - sd_score
# print(threshold)

# # Filter cells
# querry_data <- subset(querry_data, subset = prediction.score.max >= threshold)
# gc()

# DimPlot(querry_data, group.by = "class_label", label.size = 6) + theme(legend.position = "bottom",legend.text=element_text(size=4))

# labels <- unique(querry_data@meta.data$`class_label`)

# for (label in labels) {
#   plot <- DimPlot(querry_data, group.by = "class_label", cells.highlight = WhichCells(querry_data, expression = class_label == label)) + ggtitle(label)
#   print(plot)
# }

# #========================================================================
# # Step 2 : annotate with subclass
# #========================================================================

# ref_data <- readRDS(file = file.path(OUTPUTDIR, STEP_REMI, "P30_WT_remi.rds"))
# premier_modele <- ref_data@assays[["SCT"]]@SCTModel.list[[1]]
# ref_data@assays[["SCT"]]@SCTModel.list <- list(premier_modele)

# ### Set the default assay to SCT for booth seurat object
# DefaultAssay(ref_data) <- "SCT"
# DefaultAssay(querry_data) <- "SCT"

# querry_data$supertype_label <- NA

# for (class in unique(querry_data$class_label)) {
#   print("The supertype class :")
#   print(class)
  
#   # Set idents
#   Idents(querry_data) <- "class_label"
#   Idents(ref_data) <- "class_label"
  
#   # Per class
#   subset_querry_data <- subset(querry_data, idents = class)
#   subset_ref_data <- subset(ref_data, idents = class)
  
#   # Check subset sizes
#   print(paste("Number of cells in query subset:", ncol(subset_querry_data)))
#   print(paste("Number of cells in reference subset:", ncol(subset_ref_data)))

#   # Keep the same genes
#   ref_genes <- as.data.frame(rownames(subset_ref_data))
#   querry_genes <- as.data.frame(rownames(subset_querry_data))
#   all_genes <- intersect(ref_genes$`rownames(subset_ref_data)`, querry_genes$`rownames(subset_querry_data)`)
  
#   if(ncol(subset_querry_data) < K_SCORE){

#   # Calculate a new k_score
#   k_score_reference <- calculate_k_score(num_cells = ncol(subset_ref_data))
#   k_score_query <- calculate_k_score(num_cells = ncol(subset_querry_data))
#   k_score <- min(k_score_reference, k_score_query)

#   anchors <- FindTransferAnchors(reference = subset_ref_data,
#                                query = subset_querry_data,
#                                k.score = k_score,
#                                normalization.method = "SCT",
#                                features = all_genes)
# } else{
#   anchors <- FindTransferAnchors(reference = subset_ref_data,
#                                query = subset_querry_data,
#                                normalization.method = "SCT",
#                                features = all_genes)
# }

# # Number of anchors
# num_anchors <- nrow(anchors@anchors)

# if(ncol(subset_querry_data) < num_anchors){

#   # Calculate the k weight
#   k_weight <- calculate_k_weight(sample_size = num_anchors)

#   # Prediction with transfertData
#   prediction <- TransferData(anchorset = anchors,
#                            refdata = subset_ref_data$supertype_label,
#                            k.weight = k_weight
#                            )
# } else{
#   prediction <- TransferData(anchorset = anchors,
#                            refdata = subset_ref_data$supertype_label)
# }
  
#   # Add to the metadata
#   subset_querry_data <- AddMetaData(subset_querry_data, metadata = prediction)
#   querry_data$supertype_label[Cells(subset_querry_data)] <- subset_querry_data$predicted.id
# }

# ### Plot UMAP
# Idents(querry_data) <- querry_data$supertype_label
# DimPlot(querry_data, group.by = "supertype_label", label.size = 6) + theme(legend.position = "bottom",legend.text=element_text(size=4))

# ### Plot UMAP per labels
# labels <- unique(querry_data@meta.data$`supertype_label`)

# for (label in labels) {
#   plot <- DimPlot(querry_data, group.by = "supertype_label", cells.highlight = WhichCells(querry_data, expression = supertype_label == label)) + ggtitle(label)
#   print(plot)
# }

# ### UMAP plot with the prediction
# FeaturePlot(querry_data, features = "prediction.score.max", reduction = "umap")

# saveRDS(querry_data, file = "Rplots.rds")

################################################################################################ CREATE THE ARLOTTA REFERENCE SEURAT OBJECT (USED IN PREVIOUS STEPS)

# # Read the gene expression matrix (genes as columns, cells as rows)
# expression_matrix <- read.csv("/scratch/lchabot/inmed_deChevigny_scrnaseq/scRNAseq_analysis_workflow/01_Reference/reference_matrix.csv", row.names = 1)

# expression_matrix <- t(expression_matrix)

# # Read the metadata (cells as rows)
# metadata <- read.csv("/scratch/lchabot/inmed_deChevigny_scrnaseq/scRNAseq_analysis_workflow/01_Reference/metadata.csv", row.names = 1)

# cell_id_column <- "NEW_NAME"
# rownames(metadata) <- metadata[[cell_id_column]]
# metadata <- metadata[, !names(metadata) %in% cell_id_column]
# metadata <- metadata[rownames(metadata) %in% colnames(expression_matrix), ]
# expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% rownames(metadata)]

# # Create the Seurat object
# seurat_object <- CreateSeuratObject(counts = expression_matrix)

# # Add metadata to the Seurat object
# seurat_object <- AddMetaData(seurat_object, metadata = metadata)

# # Check the Seurat object
# print("Metadata in Seurat object:")
# print(head(seurat_object@meta.data))

# saveRDS(seurat_object, "/scratch/lchabot/inmed_deChevigny_scrnaseq/scRNAseq_analysis_workflow/01_Reference/modify_arlotta_seurat_object.rds")