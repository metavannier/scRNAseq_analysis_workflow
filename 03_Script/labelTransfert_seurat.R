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
library(DoubletFinder)

options(future.globals.maxSize = 3000 * 1024^2)

#=====================================
# Path to file and folder
#=====================================
DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
REF = file.path((DIRECTORY), "01_Reference")

SAMPLE_ID = "P30_WT"

ALLEN = TRUE

STEP2 = "02_seurat/"

# #========================================================================
# # Load the data an the reference
# #========================================================================

## Our data
so <- readRDS(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "_filtered_seurat_object.rds")))

if(ALLEN == TRUE){
  ref <- readRDS(file = file.path(REF, "/allen_seurat_object_with_metadata.rds"))
  ref <- PercentageFeatureSet(ref, pattern = "^mt-", col.name = "percent.mt")
  ref <- SCTransform(ref, vars.to.regress = "percent.mt", verbose = FALSE, future.seed = FALSE)

} else {
  ## The arlotta reference
  ref <- readRDS(file = file.path(REF, "/modify_arlotta_seurat_object.rds"))

  # Keep only some cells in the reference
  pattern <- "E14"
  metadata <- ref@meta.data
  filtered_metadata <- metadata[grepl(pattern, metadata$biosample_id), ]
  cells_to_keep <- rownames(filtered_metadata)
  ref <- subset(ref, cells = cells_to_keep)
  ref <- PercentageFeatureSet(ref, pattern = "^mt-", col.name = "percent.mt")
  ref <- SCTransform(ref, vars.to.regress = "percent.mt", verbose = FALSE)
}

#========================================================================
# Doublet Finder
#========================================================================

pct <- so[["pca"]]@stdev / sum(so[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

sweep.res <- paramSweep(so, PCs = 1:pcs, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
# bcmvn <- find.pK(sweep.stats)
homotypic.prop <- modelHomotypic(so$seurat_clusters)
nExp_poi <- round(0.075 * ncol(so))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
so <- doubletFinder(so, PCs = 1:pcs, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)

metadata_columns <- colnames(so@meta.data)
pattern <- "DF.classifications_0.25_0.09_\\d+"  # Adjust the regex as needed
matched_column <- metadata_columns[grepl(pattern, metadata_columns)]
num_doublets <- sub("DF\\.classifications_0\\.25_0\\.09_(\\d+)", "\\1", matched_column)
print(num_doublets)
full_column_name <- paste0("DF.classifications_0.25_0.09_", num_doublets)
cells_to_keep <- rownames(so@meta.data[so@meta.data[[full_column_name]] == "Singlet", ])
so <- subset(so, cells = cells_to_keep)

#========================================================================
# Annotation
#========================================================================

so$celltype_label <- NA
all_genes <- intersect(rownames(ref), rownames(so))
anchors <- FindTransferAnchors(reference = ref, query = so, dims = 1:30, normalization.method = "SCT", features = all_genes)
prediction <- TransferData(anchorset = anchors, refdata = as.vector(ref$celltype_label))
so <- AddMetaData(so, metadata = prediction)
so$celltype_label <- so$predicted.id

pdf(file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "_annotation_before_filtering_UMAP.pdf")), width = 8, height = 8)
plot <- DimPlot(so, group.by = "celltype_label", label.size = 6) + theme(legend.position = "bottom", legend.text = element_text(size = 4))
print(plot)

# UMAP plot per label
labels <- unique(so$celltype_label)
for (label in labels) {
  cells_to_highlight <- rownames(so@meta.data[so@meta.data$celltype_label == label, ])
  plots <- DimPlot(so, group.by = "celltype_label", cells.highlight = cells_to_highlight) +
           ggtitle(label) +
           theme(legend.position = "bottom", legend.text = element_text(size = 4))
  print(plots)
}
dev.off()

#========================================================================
# Filtering bad predicted cells : KMeans
#========================================================================

# set.seed(42)  # For reproducibility

# metadata <- so@meta.data
# metadata$cell_id <- rownames(metadata)
# unique_labels <- unique(metadata$celltype_label)
# thresholds <- data.frame(celltype_label = unique_labels, cluster_threshold = NA)

# for (label in unique_labels) {
#   subset_data <- metadata[metadata$celltype_label == label, ]
#   print(label)
#   print(nrow(subset_data))
#   # Check if there are at least 2 cells for k-means clustering
#   if (nrow(subset_data) <= 2) {
#     # If fewer than 2 cells, set the threshold to 0 (no filtering)
#     cluster_threshold <- 0
#   } else {
#     # Perform k-means clustering
#     kmeans_result <- kmeans(subset_data$prediction.score.max, centers = 2)
#     cluster_centers <- kmeans_result$centers
#     high_confidence_cluster <- which.max(cluster_centers)
#     cluster_threshold <- min(subset_data$prediction.score.max[kmeans_result$cluster == high_confidence_cluster])
#   }
  
#   thresholds[thresholds$celltype_label == label, "cluster_threshold"] <- cluster_threshold
# }

# # Merge the thresholds with the metadata
# metadata <- merge(metadata, thresholds, by = "celltype_label")

# # Apply the filtering based on the prediction score and cluster threshold
# metadata <- metadata[metadata$prediction.score.max >= metadata$cluster_threshold, ]

# saveRDS(so, file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "annotate_seurat_objet_before_filtering.rds")))

# so_clean <- subset(so, cells = metadata$cell_id)


# pdf(file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "_annotation_after_filtering_UMAP.pdf")), width = 8, height = 8)
# plot <- DimPlot(so_clean, group.by = "celltype_label", label.size = 6) + theme(legend.position = "bottom", legend.text = element_text(size = 4))
# print(plot)

# # UMAP plot per label
# labels <- unique(so_clean$celltype_label)
# for (label in labels) {
#   cells_to_highlight <- rownames(so_clean@meta.data[so_clean@meta.data$celltype_label == label, ])
#   plots <- DimPlot(so_clean, group.by = "celltype_label", cells.highlight = cells_to_highlight) +
#            ggtitle(label) +
#            theme(legend.position = "bottom", legend.text = element_text(size = 4))
#   print(plots)
# }
# dev.off()

# saveRDS(so_clean, file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "annotate_seurat_objet_after_filtering.rds")))

#========================================================================
# Filtering bad predicted cells : Remi's
#========================================================================

# 1. Get all the prediction score column names, excluding "prediction.score.max"
prediction_columns <- grep("^prediction\\.score\\.(?!max$)", colnames(so@meta.data), value = TRUE, perl = TRUE)

# 2. Create a new column "filter_celltype_label" to store the filtered labels
so@meta.data$filter_celltype_label <- "Nothing"  # Initialize with NA

# 3. Loop through each row and perform the filtering based on the criteria
for (i in 1:nrow(so@meta.data)) {
  # Extract the prediction scores for the current cell as a numeric vector
  scores <- as.numeric(so@meta.data[i, prediction_columns])

  # Sort the scores in descending order and extract the top two
  sorted_scores <- sort(scores, decreasing = TRUE)
  best_score <- sorted_scores[1]
  second_best_score <- sorted_scores[2]

  # Get the labels corresponding to the top two scores
  best_label <- prediction_columns[which.max(scores)]  # Find the label for the best score
  second_best_label <- prediction_columns[which(scores == second_best_score)[1]]  # Find the label for the second-best score

  # 4. Check if the difference between the best and second-best score is greater than 20%
  if ((best_score - second_best_score) > 0.20) {
    # Keep the label corresponding to the best score (remove "prediction.score." prefix)
    so@meta.data$filter_celltype_label[i] <- gsub("prediction\\.score\\.", "", best_label)
  } else {
    # Otherwise, set the label as NA
    so@meta.data$filter_celltype_label[i] <- "Nothing"
  }
}

saveRDS(so, file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "annotate_seurat_objet_before_filtering.rds")))

so_clean <- subset(so, subset = filter_celltype_label != "Nothing")

pdf(file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "_annotation_after_filtering_UMAP.pdf")), width = 8, height = 8)
plot <- DimPlot(so_clean, group.by = "celltype_label", label.size = 6) + theme(legend.position = "bottom", legend.text = element_text(size = 4))
print(plot)

# UMAP plot per label
labels <- unique(so_clean$celltype_label)
for (label in labels) {
  cells_to_highlight <- rownames(so_clean@meta.data[so_clean@meta.data$celltype_label == label, ])
  plots <- DimPlot(so_clean, group.by = "celltype_label", cells.highlight = cells_to_highlight) +
           ggtitle(label) +
           theme(legend.position = "bottom", legend.text = element_text(size = 4))
  print(plots)
}
dev.off()

saveRDS(so_clean, file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "annotate_seurat_objet_after_filtering.rds")))































# #========================================================================
# # Create functions
# #========================================================================

# # Function for performing label transfer with optional per-class processing
# perform_label_transfer <- function(so, ref, actual_class = "actual_class", target_label = "target_label", skip_loop = FALSE, output_dir = "default/path", pdf_file = "label_transfer_plot.pdf") {
  
#   # Initialize or clear the target_label column
#   so[[target_label]] <- NA
  
#   class_labels <- as.character(so[[actual_class]][, actual_class])

#   # Ensure the output directory exists
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir, recursive = TRUE)
#   }

#   # Create the full path for the PDF file
#   pdf_file <- file.path(output_dir, pdf_file)

#   if (!is.null(pdf_file)) {
#     pdf(pdf_file, width = 8, height = 6)  # Open PDF if specified
#   }

#   if (!skip_loop) {
#     unique_classes <- unique(class_labels)

#     for (class in unique_classes) {
#       print(class)

#       # Set idents
#       Idents(so) <- actual_class
#       Idents(ref) <- actual_class

#       # Subset per class
#       subset_querry_data <- subset(so, idents = class)
#       subset_ref_data <- subset(ref, idents = class)

#       # Keep the same genes
#       all_genes <- intersect(rownames(subset_ref_data), rownames(subset_querry_data))

#       # Directly assign actual_class to target_label
#       so[[target_label]][so[[actual_class]] == class] <- class

#       # Find anchors and Transfer Data with error handling
#       tryCatch({
#         anchors <- FindTransferAnchors(reference = subset_ref_data, query = subset_querry_data, dims = 1:30, normalization.method = "SCT", features = all_genes)
#         prediction <- TransferData(anchorset = anchors, refdata = as.vector(subset_ref_data[[target_label]]))

#         # Append target_label to each prediction column name
#         colnames(prediction) <- paste0(colnames(prediction), ".", target_label)
#         subset_querry_data <- AddMetaData(subset_querry_data, metadata = prediction)
#         new_metadata <- subset_querry_data@meta.data[, grepl(paste0(".", target_label), colnames(subset_querry_data@meta.data))]
#         so@meta.data[rownames(new_metadata), colnames(new_metadata)] <- new_metadata

#         # Assign the predicted labels
#         so[[target_label]] <- subset_querry_data[[paste0("predicted.id.", target_label)]]

#       }, error = function(e) {
#         warning(paste("Error for class:", class, "- Keeping actual_class in target_label"))
#       })
#     }
#   } else {
#     # Skip looping, transfer all cells at once
#     all_genes <- intersect(rownames(ref), rownames(so))

#     # Perform label transfer without looping
#     anchors <- FindTransferAnchors(reference = ref, query = so, dims = 1:30, normalization.method = "SCT", features = all_genes)
#     prediction <- TransferData(anchorset = anchors, refdata = as.vector(ref[[target_label]]))

#     # Append target_label to the prediction columns
#     colnames(prediction) <- paste0(colnames(prediction), ".", target_label)

#     # Add prediction metadata to the Seurat object
#     so <- AddMetaData(so, metadata = prediction)

#     # Assign the predicted labels
#     so[[target_label]] <- so[[paste0("predicted.id.", target_label)]]
#   }

#   # UMAP plot
#   plot <- DimPlot(so, group.by = target_label, label.size = 6) + 
#           theme(legend.position = "bottom", legend.text = element_text(size = 4))

#   print(plot)  # Print the main plot

#   # UMAP plot per label
#   labels <- unique(so@meta.data[[target_label]])
#   for (label in labels) {
#     cells_to_highlight <- rownames(so@meta.data[so@meta.data[[target_label]] == label, ])
#     plots <- DimPlot(so, group.by = target_label, cells.highlight = cells_to_highlight) +
#              ggtitle(label) +
#              theme(legend.position = "bottom", legend.text = element_text(size = 4))
#     print(plots)
#   }

#   if (!is.null(pdf_file)) {
#     dev.off()  # Close the PDF device
#   }

#   return(so)
# }

# # Function to calculate thresholds based on prediction scores using k-means clustering
# # Function to calculate thresholds based on prediction scores using k-means clustering
# calculate_thresholds <- function(metadata, celltype_label, prediction_column) {
#   unique_labels <- unique(metadata[[celltype_label]])
#   thresholds <- data.frame(celltype_label = unique_labels, cluster_threshold = NA, stringsAsFactors = FALSE)

#   for (label in unique_labels) {
#     subset_data <- metadata[metadata[[celltype_label]] == label & !is.na(metadata[[prediction_column]]), ]

#     # Try-catch block to handle cases with insufficient data or k-means failures
#     tryCatch({
#       if (nrow(subset_data) > 1) {
#         # Perform k-means clustering with 2 centers
#         kmeans_result <- kmeans(subset_data[[prediction_column]], centers = 2)
#         cluster_centers <- kmeans_result$centers
#         high_confidence_cluster <- which.max(cluster_centers)

#         # Calculate the cluster threshold
#         cluster_threshold <- min(subset_data[[prediction_column]][kmeans_result$cluster == high_confidence_cluster], na.rm = TRUE)
#         thresholds[thresholds$celltype_label == label, "cluster_threshold"] <- cluster_threshold
#       } else {
#         # If there are not enough rows to cluster
#         stop("Not enough data for k-means clustering")
#       }
#     }, error = function(e) {
#       # Handle error: retain original label (no filtering)
#       message(paste("Error for label:", label, "-", e$message, "- Retaining original cells (no filtering)."))
#       thresholds[thresholds$celltype_label == label, "cluster_threshold"] <- -Inf  # No filtering for this label
#     })
#   }

#   return(thresholds)
# }

# # Function to filter cells by threshold
# filter_cells_by_threshold <- function(metadata, prediction_column, thresholds, celltype_label) {
#   # Merge thresholds back to metadata
#   metadata <- merge(metadata, thresholds, by = celltype_label, by.y = "celltype_label", all.x = TRUE)

#   # Check if the merge was successful
#   if (nrow(metadata) == 0) {
#     stop("Merging thresholds into metadata resulted in zero rows. Check column names and values.")
#   }

#   # Filter cells based on prediction scores and thresholds
#   filtered_metadata <- metadata[
#     is.na(metadata[[prediction_column]]) |
#     (metadata[[prediction_column]] >= metadata$cluster_threshold),
#   ]

#   return(filtered_metadata)
# }

# # Function to create UMAP plots, save them as a PDF, and return filtered results
# create_filtered_umap <- function(so, celltype_label, output_dir = "default/path", file_name = "default_plot.pdf") {
#   set.seed(42)  # For reproducibility

#   # Ensure the output directory exists
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir, recursive = TRUE)
#   }

#   # Get metadata and create cell ID column
#   metadata <- so@meta.data
#   metadata$cell_id <- rownames(metadata)

#   # Dynamically construct the column name for prediction scores
#   prediction_column <- paste0("prediction.score.max.", celltype_label)

#   # Calculate thresholds
#   thresholds <- calculate_thresholds(metadata, celltype_label, prediction_column)
#   print(thresholds)  # Log the computed thresholds
#   write.csv(thresholds, file = file.path(output_dir, paste0(celltype_label, "_threshold.csv")), row.names = FALSE)

#   # Filter cells based on thresholds
#   filtered_metadata <- filter_cells_by_threshold(metadata, prediction_column, thresholds, celltype_label)

#   # Subset the original Seurat object using the filtered metadata
#   so <- subset(so, cells = filtered_metadata$cell_id)

#   # Create a full path for the PDF file
#   pdf_file <- file.path(output_dir, file_name)

#   # Call the plotting function with the filtered Seurat object and cell type label
#   plot_umap(so, celltype_label, pdf_file)  # Save the plots directly to the specified PDF

#   # Return the updated Seurat object and thresholds used
#   return(list(filtered_object = so, thresholds = thresholds))
# }

# # Function to plot UMAP with optional PDF saving
# plot_umap <- function(so, celltype_label, pdf_file = NULL) {
#   if (!is.null(pdf_file)) {
#     pdf(pdf_file, width = 8, height = 6)  # Open a PDF device if specified
#   }

#   # Main UMAP plot
#   plot <- DimPlot(so, group.by = celltype_label, label.size = 6) +
#           theme(legend.position = "bottom", legend.text = element_text(size = 4))
  
#   print(plot)  # Print the main plot
  
#   # Generate individual UMAP plots for each label
#   labels <- unique(so@meta.data[[celltype_label]])
#   for (label in labels) {
#     cells_to_highlight <- rownames(so@meta.data[so@meta.data[[celltype_label]] == label, ])
#     plots <- DimPlot(so, group.by = celltype_label, cells.highlight = cells_to_highlight) +
#              ggtitle(label) +
#              theme(legend.position = "bottom", legend.text = element_text(size = 4))

#     print(plots)  # Print each individual plot
#   }
  
#   if (!is.null(pdf_file)) {
#     dev.off()  # Close the PDF device if it was opened
#   }
# }

# #========================================================================
# # Load the data an the reference
# #========================================================================

# ## Our data
# so <- readRDS(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "_filtered_seurat_object.rds")))

# ## The reference
# ref <- readRDS(file = file.path(REF, "/modify_arlotta_seurat_object.rds"))

# # Keep only some cells in the reference
# pattern <- "P4"
# metadata <- ref@meta.data
# filtered_metadata <- metadata[grepl(pattern, metadata$biosample_id), ]
# cells_to_keep <- rownames(filtered_metadata)
# ref <- subset(ref, cells = cells_to_keep)
# ref <- PercentageFeatureSet(ref, pattern = "^mt-", col.name = "percent.mt")
# ref <- SCTransform(ref, vars.to.regress = "percent.mt", verbose = FALSE)

# #========================================================================
# # Doublet Finder
# #========================================================================

# pct <- so[["pca"]]@stdev / sum(so[["pca"]]@stdev) * 100
# cumu <- cumsum(pct)
# co1 <- which(cumu > 90 & pct < 5)[1]
# co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# pcs <- min(co1, co2)

# sweep.res <- paramSweep(so, PCs = 1:pcs, sct = TRUE)
# sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
# # bcmvn <- find.pK(sweep.stats)
# homotypic.prop <- modelHomotypic(so$seurat_clusters)
# nExp_poi <- round(0.075 * ncol(so))
# nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
# so <- doubletFinder(so, PCs = 1:pcs, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)

# metadata_columns <- colnames(so@meta.data)
# pattern <- "DF.classifications_0.25_0.09_\\d+"  # Adjust the regex as needed
# matched_column <- metadata_columns[grepl(pattern, metadata_columns)]
# num_doublets <- sub("DF\\.classifications_0\\.25_0\\.09_(\\d+)", "\\1", matched_column)
# print(num_doublets)
# full_column_name <- paste0("DF.classifications_0.25_0.09_", num_doublets)
# cells_to_keep <- rownames(so@meta.data[so@meta.data[[full_column_name]] == "Singlet", ])
# so <- subset(so, cells = cells_to_keep)

# #========================================================================
# # Step 1 : annotate with broad class
# #========================================================================
# so <- perform_label_transfer(so, ref, actual_class = "class_label", target_label = "class_label", skip_loop = TRUE, output_dir = file.path(OUTPUTDIR, STEP2, SAMPLE_ID), pdf_file = "class_lablel_before_filtering.pdf")
# filtered_results <- create_filtered_umap(so, celltype_label = "class_label", output_dir = file.path(OUTPUTDIR, STEP2, SAMPLE_ID), file_name = "class_lablel_after_filtering.pdf")

# #========================================================================
# # Step 2 : annotate with subclass
# #========================================================================
# so <- perform_label_transfer(so, ref, actual_class = "class_label", target_label = "subclass_label", skip_loop = FALSE, output_dir = file.path(OUTPUTDIR, STEP2, SAMPLE_ID), pdf_file = "subclass_label_before_filtering.pdf")
# filtered_results <- create_filtered_umap(so, celltype_label = "subclass_label", output_dir = file.path(OUTPUTDIR, STEP2, SAMPLE_ID), file_name = "subclass_label_after_filtering.pdf")

# #========================================================================
# # Step 3 : annotate with celltype
# #========================================================================
# so <- perform_label_transfer(so, ref, actual_class = "subclass_label", target_label = "celltype_label", skip_loop = FALSE, output_dir = file.path(OUTPUTDIR, STEP2, SAMPLE_ID), pdf_file = "celltype_label_before_filtering.pdf")
# filtered_results <- create_filtered_umap(so, celltype_label = "celltype_label", output_dir = file.path(OUTPUTDIR, STEP2, SAMPLE_ID), file_name = "celltype_label_after_filtering.pdf")

# saveRDS(so, file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID, "annotate_seurat_objet.rds")))
