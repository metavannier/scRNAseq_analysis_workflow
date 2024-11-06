#..............................................................
# merge_rds : Merging the count of samples from rds files 
#..............................................................

#-------------------------------------
# Path to files / folders
#-------------------------------------

DIRECTORY = getwd()
OUTPUTDIR = "05_Output"
SAMPLE_ID = snakemake@params[["sample_id"]]
STEP2 = "03_sims"
STEP = "04_velocity/"
TEXT_OUTPUT = snakemake@output[["merge_rds_output"]]

#-------------------------------------
# Loading library
#-------------------------------------
library(Seurat)

## Data preparation 
# Merge the rds of each samples after seurat filtration

# Read the RDS files into a list
sep_data.list <- lapply(file.path(OUTPUTDIR,STEP2, SAMPLE_ID, paste0("filtered_assigned_seurat_object.rds")), readRDS)

# Merge rds
# SO <- merge(x = sep_data.list[[1]], y = sep_data.list[2:length(sep_data.list)], merge.data = T)

# Merge the data frames in the list
merged_data <- Reduce(function(x, y) merge(x, y), sep_data.list)


## Remove/show metadata
merged_data@meta.data[grep("pANN", colnames(merged_data@meta.data), value = TRUE)] <- NULL
merged_data@meta.data[grep("RNA_snn", colnames(merged_data@meta.data), value=TRUE)] <- NULL
merged_data@meta.data[grep("diff_prob", colnames(merged_data@meta.data), value=TRUE)] <- NULL

head(merged_data@meta.data)

# Updated Jeffrey's SeuratToURD function to import the Seurat (v3 or v4) object into URD. 
seuratToURD <- function(seurat.object) {
    if (requireNamespace("Seurat", quietly = TRUE)) {
        # Create an empty URD object
        ds <- new("URD")
        
        # Copy over data
        ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
        if(!any(dim(seurat.object@assays$RNA@data) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays$RNA@data[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
        
        # Copy over metadata
        ## TO DO - grab kmeans clustering info
        get.data <- NULL
        if (.hasSlot(seurat.object, "data.info")) { 
            get.data <- as.data.frame(seurat.object@assays$RNA@data.info)
        } else if (.hasSlot(seurat.object, "meta.data")) { 
            get.data <- as.data.frame(seurat.object@meta.data) 
        }
        if(!is.null(get.data)) {
            di <- colnames(get.data)
            m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
            discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
            gi <- di[which(discrete <= 0.015)]
            ds@meta <- get.data[,m,drop=F]
            ds@group.ids <- get.data[,gi,drop=F]
        }
        
        # Copy over var.genes
        if(length(seurat.object@assays$RNA@var.features > 0)) ds@var.genes <- seurat.object@assays$RNA@var.features
        
        # Move over tSNE projection
        if (.hasSlot(seurat.object, "tsne.rot")) {
            if(!any(dim(seurat.object@tsne.rot) == 0)) {
                ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
                colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
            }
        } else if (.hasSlot(seurat.object, "reductions")) {
            if(("tsne" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$tsne) == 0)) {
                ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
                colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
            }
        }
        
        # Move over PCA results
        if (.hasSlot(seurat.object, "pca.x")) {
            if(!any(dim(seurat.object@pca.x) == 0)) {
                ds@pca.load <- seurat.object@pca.x
                ds@pca.scores <- seurat.object@pca.rot
                warning("Need to set which PCs are significant in @pca.sig")
            }
            ## TO DO: Convert SVD to sdev
        } else if (.hasSlot(seurat.object, "reductions")) {
            if(("pca" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$pca@gene.loadings) == 0)) {
                ds@pca.load <- as.data.frame(seurat.object@reductions$pca@gene.loadings)
                ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
                ds@pca.sdev <- seurat.object@reductions$pca@sdev
                ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
            }
        }
        return(ds)
    } else {
        stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
    }
}

# Format of an URD object (because trouble when charge URD environment with Seurat environment)

URD <- methods::setClass("URD", slots = c(
    count.data=c("dgCMatrix", NULL), 
    logupx.data=c("dgCMatrix", NULL), 
    meta="data.frame", 
    group.ids="data.frame", 
    var.genes="vector", 
    knn="list",
    pca.sdev="vector", 
    pca.load="data.frame", 
    pca.scores="data.frame", 
    pca.sig="vector", 
    tsne.y="data.frame", 
    plot.3d="list",
    gene.sig.z="data.frame", 
    dm=c("DiffusionMap",NULL), 
    diff.data="data.frame", 
    pseudotime="data.frame",
    pseudotime.stability="list", 
    tree="list",
    nmf.g=c("dgCMatrix", "matrix", NULL),
    nmf.c=c("dgCMatrix", "matrix", NULL)
))

urdO <- seuratToURD(merged_data)
saveRDS(urdO, file = file.path( OUTPUTDIR, STEP, "merged_urd_object.rds"))


saveRDS(merged_data, file = file.path( OUTPUTDIR, STEP, "merged_seurat_object.rds"))

# Extract raw count data
count_data <- GetAssayData(merged_data, slot = "counts")
# Save count data to a text file
write.table(as.matrix(count_data), file = file.path(OUTPUTDIR,STEP, paste0("merge_count_data.txt")), sep = "\t", 
            quote = FALSE, col.names = NA)  # Use col.names = NA to not include column names

# Extract metadata
metadata <- merged_data@meta.data
# Save metadata to a text file
write.table(metadata, file = file.path(OUTPUTDIR,STEP, paste0("merge_metadata.txt")), sep = "\t", 
            quote = FALSE, row.names = TRUE)  # Set row.names = TRUE to include row names

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

output_file<-file(TEXT_OUTPUT)
writeLines(c("Merging step finished"), output_file)
close(output_file)

sessionInfo(package = NULL)