#..............................................................
# Preparation of the folowing reference for sims :
# Paola Arlotta reference : 
#..............................................................

#-------------------------------------
# Library
#-------------------------------------
library(rhdf5)
library(data.table)
library(Seurat)
library(dplyr)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
REF = file.path((DIRECTORY), "01_Reference")

SAMPLE_ID = snakemake@params[["sample_id"]]
REFERENCE_NAME = snakemake@params[["reference_name"]]

NORM_METHOD = snakemake@params[["norm_method"]]
NORM_SCALE_FACTOR = snakemake@params[["norm_scale_factor"]]

ALLEN_METADATA = snakemake@params[["allen_metadata"]]
ALLEN_MATRIX = snakemake@params[["allen_matrix"]]

OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]
 
STEP2 = "02_seurat/"
STEP3 = "03_sims/"
STEP_REMI = "03_remi/"

METADATA_REMI <- TRUE

#-------------------------------------
# Library
#-------------------------------------
library(rhdf5)
library(data.table)
library(Seurat)
library(dplyr)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
REF = file.path((DIRECTORY), "01_Reference")

SAMPLE_ID = snakemake@params[["sample_id"]]
REFERENCE_NAME = snakemake@params[["reference_name"]]

NORM_METHOD = snakemake@params[["norm_method"]]
NORM_SCALE_FACTOR = snakemake@params[["norm_scale_factor"]]

ALLEN_METADATA = snakemake@params[["allen_metadata"]]
ALLEN_MATRIX = snakemake@params[["allen_matrix"]]

OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]
 
STEP2 = "02_seurat/"
STEP3 = "03_sims/"
STEP_REMI = "03_remi/"

METADATA_REMI <- TRUE


metadata <- fread(file = file.path(REF, ALLEN_METADATA))
print("Reference metadata loaded. It's size is :")
print(dim(metadata))

metadata <- metadata[metadata$region_label %in% "SSp",]
print("Only SSp are kept in the metadata: The new size of this one is :")
print(dim(metadata))

# Load the HDF5 matrix
h5_file <- H5Fopen(file.path(REF, ALLEN_MATRIX))
h5ls(h5_file)
expression_matrix <- h5read(h5_file, "/data/counts")
gene_names <- h5read(h5_file, "/data/gene")
sample_names <- h5read(h5_file, "/data/samples")
H5Fclose(h5_file)

colnames(expression_matrix) <- as.character(gene_names)  # Si nécessaire
rownames(expression_matrix) <- as.character(sample_names)  # Si nécessaire

expression_matrix <- expression_matrix[rownames(expression_matrix) %in% metadata$sample_name,]
metadata <- metadata[metadata$sample_name %in% rownames(expression_matrix),]

print("The metadata and matrix of reference have the same cells now : The new size of them two :")
print(dim(expression_matrix))
print(dim(metadata))

gc()

expression_matrix <- t(expression_matrix)
print("The matrix is transposed")

gc()

colnames(expression_matrix) <- gsub("_", "-", original_sample_names)

seurat_object <- CreateSeuratObject(counts = expression_matrix, project = "allen_seuratObject")
saveRDS(seurat_object, file = file.path(REF, "allen_seurat_object.rds"))

rm("seurat_object")

gc()

#### Pour une raison que j'ignore cette partie ne marche pas, je dois le faire en local

seurat_object <- readRDS(file = file.path(REF, "allen_seurat_object.rds"))

metadata$sample_name <- gsub("_", "-", metadata$sample_name)
metadata <- metadata[match(colnames(seurat_object), metadata$sample_name), ]
seurat_object <- AddMetaData(object = seurat_object, metadata = metadata)

saveRDS(seurat_object, file = file.path(REF, "allen_seurat_object_with_metadata.rds"))



# #-------------------------------------
# # Load the metadata of the reference
# #-------------------------------------

# if(METADATA_REMI){
#     reference_metadata <- fread(file = file.path(OUTPUTDIR, STEP_REMI, "final_allen.csv"))
#     cells_names <- "sample_name"
#     celltypes <- "celltype_label"
# } else {
#     reference_metadata <- fread(file = file.path(REF, ALLEN_METADATA))
#     cells_names <- "sample_name"
#     celltypes <- "subclass_label"
# }

# print("Reference metadata loaded. It's size is :")
# print(dim(reference_metadata))

# ########################################## PROPRE AU PROJET ##########################################
# #-------------------------------------
# # Keep only some part of the cortex
# # we want to study: Ssp
# #-------------------------------------
# if(!METADATA_REMI){
#     reference_metadata <- reference_metadata[reference_metadata$region_label %in% "SSp",]
#     print("Only SSp are kept in the metadata: The new size of this one is :")
#     print(dim(reference_metadata))
# }
# ######################################################################################################

# #-------------------------------------
# # For the training we need at least
# # two cells per labels
# #-------------------------------------
# occurence <- table(reference_metadata[[celltypes]])
# one_cell_per_label <- names(occurence[occurence == 1])
# reference_metadata <- reference_metadata[!(reference_metadata[[celltypes]] %in% one_cell_per_label), ]

# print("The new size after removing labels thaht only appear once :")
# print(dim(reference_metadata))

# #-------------------------------------
# # Load the matrix of reference
# #-------------------------------------
# reference_matrix <- h5read(file = file.path(REF, ALLEN_MATRIX), "/data/counts")
# reference_genes <- h5read(file = file.path(REF, ALLEN_MATRIX), "/data/gene")
# reference_cells <- h5read(file = file.path(REF, ALLEN_MATRIX), "/data/samples")

# rownames(reference_matrix) <- as.character(reference_cells)
# colnames(reference_matrix) <- as.character(reference_genes)

# rm(list = c("reference_genes","reference_cells"))
# gc()

# print("Reference matrix loaded. The size is :")
# print(dim(reference_matrix))

# #-------------------------------------
# # Keep cells that match both the 
# # metadata and matrix of reference
# #-------------------------------------
# reference_matrix <- reference_matrix[rownames(reference_matrix) %in% reference_metadata[[cells_names]],]
# reference_metadata <- reference_metadata[reference_metadata[[cells_names]] %in% rownames(reference_matrix),]

# print("The metadata and matrix of reference have the same cells now : The new size of them two :")
# print(dim(reference_matrix))
# print(dim(reference_metadata))
# print(reference_matrix[1:5, 1:25])

# #-------------------------------------
# # Order metadata according to matrix
# # (Necessary ?)
# #-------------------------------------
# reference_metadata <- reference_metadata[match(rownames(reference_matrix),reference_metadata[[cells_names]]),]
# print("Cells are in the same order in both files")

# #-------------------------------------
# # Write the the metadata as csv   
# #-------------------------------------
# dir.create(file.path(OUTPUTDIR, STEP3, REFERENCE_NAME))

# fwrite(x = reference_metadata, file = file.path(OUTPUTDIR, STEP3, REFERENCE_NAME, paste0(OUTPUT_NAME_REF_METADATA, ".csv")))

# rm("reference_metadata")
# gc()

# print("Metadata is written")

# #-------------------------------------
# # Normalization of the matrix of
# # reference
# #-------------------------------------
# reference_matrix <- NormalizeData( object = reference_matrix, normalization.method = NORM_METHOD, scale.factor = NORM_SCALE_FACTOR, verbose = FALSE)
# reference_matrix <- as.matrix(reference_matrix)

# print("The matrix of reference is now normalized the same way as our matrix to analyze: Here the first 10 rows and 10 columns :")
# print(reference_matrix[1:5, 1:25])

# #-------------------------------------
# # Load our matrix to annotate
# #-------------------------------------
# dir.create(file.path(OUTPUTDIR, STEP3, SAMPLE_ID))

# matrix <- fread(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_normalized_matrix.csv")))
# matrix <- as.matrix(matrix)
# rownames(matrix) <- matrix[,"V1"]
# matrix <- matrix[, -1]

# print("Our matrix to annotate is loaded, see the dimension below :")
# print(dim(matrix))

# # -------------------------------------
# # Wtrite the two new matrix
# # -------------------------------------
# class(matrix) <- "numeric"
# write.csv(matrix, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_MATRIX,".csv")))
# print("Our matrix is written")

# rm("matrix")
# gc()

# class(reference_matrix) <- "numeric"
# write.csv(reference_matrix, file.path(OUTPUTDIR, STEP3, REFERENCE_NAME, paste0(OUTPUT_NAME_REF_MATRIX,".csv")))
# print("Our matrix of reference is written")

# rm("reference_matrix")
# gc()

#-------------------------------------
# create the output file
#-------------------------------------
output_file<-file(TEXT_OUTPUT)
writeLines(c("Files preparation for the Allen reference for SIMS finished (CSV format)"), output_file)
close(output_file)





# #-------------------------------------
# # Look at genes who are common in both
# # matrix and subset genes who are nots
# #-------------------------------------
# gene_intersection <- intersect(colnames(matrix), colnames(reference_matrix))
# reference_matrix <- reference_matrix[, gene_intersection]
# matrix <- matrix[, gene_intersection]

# print("Only the gene in common in booth our matrix are kept : The new dimmension :")
# print(dim(reference_matrix))
# print(colnames(reference_matrix)[1:35])
# print(dim(matrix))
# print(colnames(matrix)[1:35])