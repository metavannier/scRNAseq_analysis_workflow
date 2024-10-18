# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

# ..........................................................................................................
## @knitr HTODemux
# ..........................................................................................................

# if(as.logical(MULTIPLEX) == TRUE){
#   # Demultiplexing data 
#   seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = 0.99)
#   cat("<br><br>Removed cells after filtering Doublet et Negative:", sum(seurat_obj[["hash.ID"]] == "Negative" | seurat_obj[["hash.ID"]] == "Doublet"));
#   cat("<br><br>Remaining cells after filtering:", sum(seurat_obj$hash.ID %in% HTO));
#   cat("\n<br>\n");
#   seurat_obj <- subset(seurat_obj, idents = c("Doublet","Negative"), invert = TRUE) ### Keep only the Singlet cells (Remove Negative and Doublet)

#   seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "mean.var.plot")
#   seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
#   seurat_obj <- RunPCA(seurat_obj)
#   seurat_obj <- RunUMAP(seurat_obj, dims = 1:25)
#   print(DimPlot(seurat_obj))
# }

# ..........................................................................................................
## @knitr ldr
# ..........................................................................................................

# Perform linear dimensional reduction
seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = FALSE)
VizDimLoadings(seurat_obj, dims = 1:DIMS, reduction = "pca")
# DimPlot(seurat_obj, reduction = "pca", group.by = "orig.ident")
# DimPlot(seurat_obj, reduction = "pca", group.by = "categorie")
DimHeatmap(seurat_obj, dims = 1:DIMS, cells = 500, balanced = TRUE)

# ..........................................................................................................
## @knitr elbowplot
# ..........................................................................................................

ElbowPlot(seurat_obj, ndims = DIMS)

# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))

# Elbow plot to visualize 
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
  
# ..........................................................................................................
## @knitr heterogeneity_pca
# ..........................................................................................................

# Compute PCA on selected variable genes
seurat_obj <- RunPCA( object = seurat_obj,
                 npcs     = pcs,
                 verbose  = .VERBOSE,
                 assay = "SCT");
