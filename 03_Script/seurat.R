## @knitr qcmetrics

# Load output
OUTPUT = snakemake@output[["seurat_object"]]
COUNT = snakemake@output[["count_matrix"]]
DATA = snakemake@output[["data_matrix"]]
SCALE = snakemake@output[["scale_data_matrix"]]
# Load parameters
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
EXPCLUST = snakemake@params[["expclust"]]
CLUSTERIDS = snakemake@params[["cluster_ids"]]

# Load the single cell dataset
sc_data <- Read10X(data.dir = SC_DATA_PATH)
# Initialize the Seurat object with the raw (non-normalized data).
obj_seurat <- CreateSeuratObject(counts = sc_data, project = SAMPLE, min.cells = MINCELLS, min.features = MINFEATURES)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
obj_seurat[["percent.mt"]] <- PercentageFeatureSet(obj_seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(obj_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter to visualize mt-RNA and feature-RNA relationships
plot1 <- FeatureScatter(obj_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## Cell filtration
# Filter genes
# The gene should be expressed at least 1% of samples
obj_seurat <- subset(obj_seurat, subset = nFeature_RNA > MINFEATURERNA & nFeature_RNA < MAXFEATURERNA & percent.mt < PERCENTMT)

## @knitr hvf

### Workflow for data transformation
if (NORMWF == "basictransform"){
    ## Normalize the data
    obj_seurat <- NormalizeData(obj_seurat, normalization.method = NORMMETHOD, scale.factor = SCALEFACTOR)

    ## Identification of highly variable features
    obj_seurat <- FindVariableFeatures(obj_seurat, selection.method = SELECTMETHOD, nfeatures = NFEATURES)

    # Identify the 10 most highly variable genes
    top <- head(VariableFeatures(obj_seurat), NHVG)

        # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(obj_seurat)
    plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)
    plot2

        # Scaling the data and remove unwanted sources of variation (mitochondrial)
    all.genes <- rownames(obj_seurat)
    obj_seurat <- ScaleData(obj_seurat, features = all.genes, vars.to.regress = "percent.mt")
    # USE ONLY FOR THE TEST (select the nb of variable genes for the scale) BUT NOT USE FOR THE HEATMAP
    #obj_seurat <- ScaleData(obj_seurat)
} else if (NORMWF == "SCTransform"){        
        # run sctransform
        obj_seurat <- SCTransform(obj_seurat, new.assay.name = "SCT", method = "glmGamPoi", variable.features.n = NFEATURES, vars.to.regress = "percent.mt", verbose = FALSE)
        # Identify the 10 most highly variable genes
        top <- head(VariableFeatures(obj_seurat), NHVG)
        plot1 <- VariableFeaturePlot(obj_seurat)
        plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)
        plot2
} else {
        print("You have to enter the normalization workflow to use in the config.yaml file. Could be SCTransform or basictransform")
}

## @knitr ldr

# Perform linear dimensional reduction
obj_seurat <- RunPCA(obj_seurat, features = VariableFeatures(object = obj_seurat))
#obj_seurat <- RunPCA(obj_seurat, verbose = FALSE)
VizDimLoadings(obj_seurat, dims = 1:10, reduction = "pca")

DimPlot(obj_seurat, reduction = "pca")

DimHeatmap(obj_seurat, dims = 1:DIMS, cells = 500, balanced = TRUE)

## @knitr elbowplot
ElbowPlot(obj_seurat)

## @knitr cluster
obj_seurat <- FindNeighbors(obj_seurat, dims = 1:DIMS, verbose = FALSE)
obj_seurat <- FindClusters(obj_seurat, resolution = RESOLUTION, verbose = FALSE)

# Look at cluster IDs and rename it according to your marker genes
#WhichCells(obj_seurat, idents = "0")
cluster.ids = as.list(strsplit(CLUSTERIDS, ",")[[1]])
names(cluster.ids) <- levels(obj_seurat)
obj_seurat <- RenameIdents(obj_seurat, cluster.ids)

## @knitr umap
obj_seurat <- RunUMAP(obj_seurat, dims = 1:DIMS)
DimPlot(obj_seurat, reduction = "umap")

obj_seurat <- RunTSNE(obj_seurat, dims = 1:DIMS)
DimPlot(obj_seurat, reduction = "tsne")

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
