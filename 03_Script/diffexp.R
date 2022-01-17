## @knitr violin

# Load input
SEURAT_OBJECT = snakemake@input[["seurat_object"]]
#Load path
PATH <- getwd()
PATH <- gsub("/03_Script","", PATH)
# Load output
VIOLINPLOT = snakemake@output[["violinplot"]]
UMAPFEATURE = snakemake@output[["umapfeature"]]
TSNEFEATURE = snakemake@output[["tsnefeature"]]
RIDGEFEATURE = snakemake@output[["ridgefeature"]]
HEATMAPFEATURE = snakemake@output[["heatmapfeature"]]
# Load parameters
FEATURES = snakemake@params[["features"]]
# Load seurat object generated after clustering
obj_seurat <- readRDS(SEURAT_OBJECT)

featureslist = c(strsplit(FEATURES, ",")[[1]])


ViolinPlot <- VlnPlot(obj_seurat, features = featureslist)
plot(ViolinPlot)

for (i in 1:length(featureslist)) {
    pdffile = paste(PATH,"/05_Output/03_diffexp/violin_plot/",featureslist[i],"_violin_plot.pdf",sep = "")
    pdf(pdffile)
    vln = VlnPlot(obj_seurat, features = featureslist[i])
    print(vln)
    dev.off()
}

## @knitr umapfeature

umapfeature <- FeaturePlot(obj_seurat, features = featureslist, reduction = "umap")
plot(umapfeature)

for (i in 1:length(featureslist)) {
    pdffile = paste(PATH,"/05_Output/03_diffexp/umap_plot/",featureslist[i],"_umapfeature_plot.pdf",sep = "")
    pdf(pdffile)
    umapfeature <- FeaturePlot(obj_seurat, features = featureslist[i], reduction = "umap")
    print(umapfeature)
    dev.off()
}

## @knitr tsnefeature

tsnefeature <- FeaturePlot(obj_seurat, features = featureslist, reduction = "tsne")
plot(tsnefeature)

for (i in 1:length(featureslist)) {
    pdffile = paste(PATH,"/05_Output/03_diffexp/tsne_plot/",featureslist[i],"_tsnefeature_plot.pdf",sep = "")
    pdf(pdffile)
    tsnefeature <- FeaturePlot(obj_seurat, features = featureslist[i], reduction = "tsne")
    print(tsnefeature)
    dev.off()
}

## @knitr ridgefeature

ridgefeature <- RidgePlot(obj_seurat, features = featureslist)
plot(ridgefeature)

for (i in 1:length(featureslist)) {
    pdffile = paste(PATH,"/05_Output/03_diffexp/ridge_plot/",featureslist[i],"_ridgefeature_plot.pdf",sep = "")
    pdf(pdffile)
    ridgefeature <- RidgePlot(obj_seurat, features = featureslist[i])
    print(ridgefeature)
    dev.off()
}

## @knitr heatmapfeature

heatmapfeature <- DoHeatmap(object = obj_seurat, features = featureslist)
plot(heatmapfeature)
pdf(HEATMAPFEATURE)
print(heatmapfeature)
dev.off()

## Pour representer volcanoplot faire :Please change the threshold to 0.0. This will mean, however, that FindMarkers() takes longer to complete.
#puis générer le volcano.
# dans FindMarkers faire comparaison de normal melanocyte contre tout les autres cell type? ou faire pour normal contre chaque type? les deux?


# ## @knitr info
sessionInfo()

