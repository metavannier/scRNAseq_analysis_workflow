## @knitr qcmetrics

# Load output
OUTPUT = snakemake@params[["seurat_object"]]
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
VizDimLoadings(obj_seurat, dims = 1:5, reduction = "pca")

DimPlot(obj_seurat, reduction = "pca")

DimHeatmap(obj_seurat, dims = 1:15, cells = 500, balanced = TRUE)

## @knitr elbowplot
ElbowPlot(obj_seurat)

## @knitr cluster
obj_seurat <- FindNeighbors(obj_seurat, dims = 1:20, verbose = FALSE)
obj_seurat <- FindClusters(obj_seurat, resolution = RESOLUTION, verbose = FALSE)

## @knitr umap
obj_seurat <- RunUMAP(obj_seurat, dims = 1:20)
DimPlot(obj_seurat, reduction = "umap")

#saveRDS(obj_seurat, file = OUTPUT)

## @knitr info
sessionInfo()

