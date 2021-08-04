## @knitr qcmetrics

# Load the single cell dataset
sc-data <- Read10X(data.dir = SC_DATA_PATH)
# Initialize the Seurat object with the raw (non-normalized data).
sc <- CreateSeuratObject(counts = sc-data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)