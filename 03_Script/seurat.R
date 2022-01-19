## @knitr qcmetrics
# Load input
METADATA = snakemake@input[["aggrcsv"]]
# Load output
OUTPUT = snakemake@output[["seurat_object"]]
COUNT = snakemake@output[["count_matrix"]]
DATA = snakemake@output[["data_matrix"]]
SCALE = snakemake@output[["scale_data_matrix"]]
# Load parameters
SAMPLES = snakemake@params[["samples"]]
WT = snakemake@params[["wt"]]
MINCELLS = snakemake@params[["min_cells"]]
MINFEATURES = snakemake@params[["min_features"]]
MINFEATURERNA = snakemake@params[["minFeature_RNA"]]
MAXFEATURERNA = snakemake@params[["maxFeature_RNA"]]
PERCENTMT = snakemake@params[["percent_mt"]]
NORMWF = snakemake@params[["normwf"]]
NORMMETHOD = snakemake@params[["normalization_method"]]
SCALEFACTOR = snakemake@params[["scale_factor"]]
NFEATURES = snakemake@params[["nfeatures"]]
SELECTMETHOD = snakemake@params[["selection_method"]]
NHVG = snakemake@params[["nHVG"]]
DIMS = snakemake@params[["dims"]]
RESOLUTION = snakemake@params[["resolution"]]
CLUSTERIDS = snakemake@params[["cluster_ids"]]

# Load the single cell dataset
sc_data <- Read10X(data.dir = SC_DATA_PATH)

# Initialize the Seurat object with the raw (non-normalized data).
obj_seurat <- CreateSeuratObject(counts = sc_data, project = NPROJ, min.cells = MINCELLS, min.features = MINFEATURES, names.field = 2, names.delim = "\\-")

# Rename the sample name of the aggregated Seurat object 
sample_list = strsplit(SAMPLES, ",")[[1]]
samplename = obj_seurat@meta.data$orig.ident

nameid = rep(sample_list[1],length(samplename))
for (i in 2:length(sample_list)) {
    nameid[samplename %in% c(i)] = c(sample_list[i])
}
names(nameid) = rownames(obj_seurat@meta.data)

obj_seurat$original <- obj_seurat$orig.ident
obj_seurat$orig.ident <- NULL
obj_seurat <- AddMetaData(
  object = obj_seurat,
  metadata = nameid,
  col.name = "orig.ident")

# Add a categorie column for WT and "mutant"
samplename = obj_seurat@meta.data$orig.ident

categorieid = rep("Melanocytes",length(samplename))
categorieid[samplename %in% c(WT)] = "Normal"
names(categorieid) = rownames(obj_seurat@meta.data)

obj_seurat <- AddMetaData(
  object = obj_seurat,
  metadata = categorieid,
  col.name = "categorie")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
obj_seurat[["percent.mt"]] <- PercentageFeatureSet(obj_seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(obj_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
# FeatureScatter to visualize mt-RNA and feature-RNA relationships
plot1 <- FeatureScatter(obj_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot2 <- FeatureScatter(obj_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot1 + plot2

## Cell filtration
# Filter genes
# The gene should be expressed at least 1% of samples
obj_seurat <- subset(obj_seurat, subset = nFeature_RNA > MINFEATURERNA & nFeature_RNA < MAXFEATURERNA & percent.mt < PERCENTMT)

## @knitr hvf

### Workflow for data transformation
### Need to include the integration procedure if cell ranger aggregate is used
if (NORMWF == "basictransform"){
    ## Normalize the data
    obj_seurat <- NormalizeData(obj_seurat, normalization.method = NORMMETHOD, scale.factor = SCALEFACTOR)

    ## Identification of highly variable features
    obj_seurat <- FindVariableFeatures(obj_seurat, selection.method = SELECTMETHOD, nfeatures = NFEATURES)

    # Identify the most highly variable genes
    top <- head(VariableFeatures(obj_seurat), NHVG)

    # Plot variable features with and without labels
    plot1 <- VariableFeaturePlot(obj_seurat)
    plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)
    plot2

    # Scaling the data and remove unwanted sources of variation (mitochondrial)
    all.genes <- rownames(obj_seurat)
    obj_seurat <- ScaleData(obj_seurat, features = all.genes, vars.to.regress = "percent.mt")
    # USE ONLY FOR THE TEST (select the nb of variable genes for the scale) BUT NOT USE FOR THE HEATMAP
    #obj_seurat <- ScaleData(obj_seurat)
    ## Normalize and correct the batch effect using the Seurat integration procedure
} else if (NORMWF == "SCTransform"){
    obj_seurat <- SCTransform(obj_seurat, new.assay.name = "SCT", method = "glmGamPoi", variable.features.n = NFEATURES, vars.to.regress = "percent.mt", verbose = FALSE)

    ### If you want to do the integration to correct the batch effect :
    # Normalize datasets individually by SCTransform()
    # obj_seurat.list <- SplitObject(obj_seurat, split.by = "orig.ident")
    # obj_seurat.list <- lapply(X = obj_seurat.list, FUN = function(x) {
    #   x <- SCTransform(object = x, new.assay.name = "SCT", method = "glmGamPoi", variable.features.n = NFEATURES, vars.to.regress = "percent.mt", verbose = FALSE)
    # })
    # features <- SelectIntegrationFeatures(object.list = obj_seurat.list, nfeatures = NFEATURES)
    # # Run the PrepSCTIntegration() function prior to identifying anchors
    # obj_seurat.list <- PrepSCTIntegration(object.list = obj_seurat.list, anchor.features = NFEATURES)
    # # Find a set of anchors between a list of Seurat objects
    # obj_seurat.anchors <- FindIntegrationAnchors(object.list = obj_seurat.list, normalization.method = "SCT", anchor.features = features)
    # # Perform dataset integration using a pre-computed AnchorSet
    # obj_seurat <- IntegrateData(anchorset = obj_seurat.anchors, normalization.method = "SCT")
    # DefaultAssay(obj_seurat) <- "integrated"

    # Identify the 10 most highly variable genes
    ## Don't work after integration
    top <- head(VariableFeatures(obj_seurat), NHVG)
    plot1 <- VariableFeaturePlot(obj_seurat)
    plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)
    plot2
} else {
    print("You have to enter the normalization workflow to use in the config.yaml file. Could be SCTransform or basictransform")
}

### Add subpopulation identity as metadata for the Delfini project ###
    ## For the project of Delphini we have to find the different melanoma cell line from the paper Wouters et al., 2019 
    # List of different object corresponding to the dataset :
    obj_seurat_list <- SplitObject(obj_seurat, split.by = "orig.ident")
    # Perform linear dimensional reduction
    obj_seurat_list$Mix_MM_lines <- RunPCA(obj_seurat_list$Mix_MM_lines, features = VariableFeatures(object = obj_seurat_list$Mix_MM_lines))
    obj_seurat_list$Mix_MM_lines <- FindNeighbors(obj_seurat_list$Mix_MM_lines, dims = 1:DIMS, verbose = FALSE)
    obj_seurat_list$Mix_MM_lines <- FindClusters(obj_seurat_list$Mix_MM_lines, resolution = RESOLUTION, verbose = FALSE)

    # Look at cluster IDs and rename it according to your marker genes
    obj_seurat_list$Mix_MM_lines <- RenameIdents(obj_seurat_list$Mix_MM_lines, `0` = "MM087", `1` = "A375",`2` = "MM029",`3` = "MM099",`4` = "MM031",`5` = "MM057",`6` = "MM001",`7` = "MM074",`8` = "MM047",`9` = "MM011")
    # colnames(x = obj_seurat_list$Mix_MM_lines)
    obj_seurat_list$Mix_MM_lines[["cluster"]] <- Idents(object = obj_seurat_list$Mix_MM_lines)
    obj_seurat_list$M12 <- RenameIdents(obj_seurat_list$M12, `1` = "M12")
    obj_seurat_list$M12[["cluster"]] <- Idents(object = obj_seurat_list$M12)
    obj_seurat_list$M15 <- RenameIdents(obj_seurat_list$M15, `2` = "M15")
    obj_seurat_list$M15[["cluster"]] <- Idents(object = obj_seurat_list$M15)
    obj_seurat_list$M27 <- RenameIdents(obj_seurat_list$M27, `3` = "M27")
    obj_seurat_list$M27[["cluster"]] <- Idents(object = obj_seurat_list$M27)
    obj_seurat_list$Normal_Melanocytes <- RenameIdents(obj_seurat_list$Normal_Melanocytes, `4` = "Normal_Melanocytes")
    obj_seurat_list$Normal_Melanocytes[["cluster"]] <- Idents(object = obj_seurat_list$Normal_Melanocytes)

    obj_seurat <- merge(obj_seurat_list$M12,y = c(obj_seurat_list$M15,obj_seurat_list$M27,obj_seurat_list$Normal_Melanocytes,obj_seurat_list$Mix_MM_lines))
    obj_seurat <- SCTransform(obj_seurat, new.assay.name = "SCT", method = "glmGamPoi", variable.features.n = NFEATURES, vars.to.regress = "percent.mt", verbose = FALSE)
    ### End of the code for Delfini project

## @knitr ldr
# Perform linear dimensional reduction
obj_seurat <- RunPCA(obj_seurat, features = VariableFeatures(object = obj_seurat))
VizDimLoadings(obj_seurat, dims = 1:10, reduction = "pca")
DimPlot(obj_seurat, reduction = "pca", group.by = "orig.ident")
DimPlot(obj_seurat, reduction = "pca", group.by = "categorie")
DimHeatmap(obj_seurat, dims = 1:DIMS, cells = 500, balanced = TRUE)

## @knitr elbowplot
ElbowPlot(obj_seurat)

## @knitr cluster
### We already have the cluster for this analyse
obj_seurat <- RunPCA(obj_seurat, features = VariableFeatures(object = obj_seurat))
# obj_seurat <- FindNeighbors(obj_seurat, dims = 1:DIMS, verbose = FALSE)
# obj_seurat <- FindClusters(obj_seurat, resolution = RESOLUTION, verbose = FALSE)

# # Look at cluster IDs and rename it according to your marker genes
# cluster.ids = as.list(strsplit(CLUSTERIDS, ",")[[1]])
# names(cluster.ids) <- levels(obj_seurat)
# obj_seurat <- RenameIdents(obj_seurat, cluster.ids)

## @knitr umap
obj_seurat <- RunUMAP(obj_seurat, dims = 1:DIMS)
dittoDimPlot(obj_seurat, reduction = "umap", var = "cluster", do.label = TRUE)
DimPlot(obj_seurat, reduction = "umap", group.by = "categorie")
obj_seurat <- RunTSNE(obj_seurat, dims = 1:DIMS)
dittoDimPlot(obj_seurat, reduction = "tsne", var = "cluster", do.label = TRUE)
DimPlot(obj_seurat, reduction = "tsne", group.by = "categorie")

## @knitr datasave
# Saving the seurat object
saveRDS(obj_seurat, file = OUTPUT)

# How many cells are in each cluster
print("\nNumber of cells in each cluster\n")
table(Idents(obj_seurat))

# Retrieve data in a count matrix (counts)
counts <- GetAssayData(object = obj_seurat, slot = "counts")
counts <- as(Class = 'matrix', object = counts)
write.csv(x = counts, file = COUNT, quote = FALSE)

# Retrieve data in a nonnormalized expression matrix (data)
data <- GetAssayData(object = obj_seurat)
data <- as(Class = 'matrix', object = data)
write.csv(x = data, file = DATA, quote = FALSE)

# Retrieve data in a normalized expression matrix (scale)
scale <- GetAssayData(object = obj_seurat, slot = "scale.data")
scale <- as(Class = 'matrix', object = scale)
write.csv(x = scale, file = SCALE, quote = FALSE)

## @knitr info
session_info()
