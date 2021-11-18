## @knitr violin

# Load input
SEURAT_OBJECT = snakemake@input[["seurat_object"]]
# Load output
VIOLINPLOT = snakemake@output[["violinplot"]]
UMAPFEATURE = snakemake@output[["umapfeature"]]
TSNEFEATURE = snakemake@output[["tsnefeature"]]
RIDGEFEATURE = snakemake@output[["ridgefeature"]]
HEATMAPFEATURE = snakemake@output[["heatmapfeature"]]
# Load parameters
FEATURES = snakemake@params[["features"]]
# Load seurat object generated after clustering
seurat_object <- readRDS(SEURAT_OBJECT)

featureslist = c(strsplit(FEATURES, ",")[[1]])
ViolinPlot <- VlnPlot(seurat_object, features = featureslist)
plot(ViolinPlot)
pdf(VIOLINPLOT)
print(ViolinPlot)
dev.off()

## @knitr umapfeature

umapfeature <- FeaturePlot(seurat_object, features = featureslist, reduction = "umap")
plot(umapfeature)
pdf(UMAPFEATURE)
print(umapfeature)
dev.off()

## @knitr tsnefeature

tsnefeature <- FeaturePlot(seurat_object, features = featureslist, reduction = "tsne")
plot(tsnefeature)
pdf(TSNEFEATURE)
print(tsnefeature)
dev.off()

## @knitr ridgefeature

ridgefeature <- RidgePlot(seurat_object, features = featureslist)
plot(ridgefeature)
pdf(RIDGEFEATURE)
print(ridgefeature)
dev.off()

## @knitr heatmapfeature

heatmapfeature <- DoHeatmap(object = seurat_object, features = featureslist)
plot(heatmapfeature)
pdf(HEATMAPFEATURE)
print(heatmapfeature)
dev.off()

# ## @knitr info
# sessionInfo()

