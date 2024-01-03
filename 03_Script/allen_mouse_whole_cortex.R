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

NORM_METHOD = snakemake@params[["norm_method"]]
NORM_SCALE_FACTOR = snakemake@params[["norm_scale_factor"]]

ALLEN_METADATA = snakemake@params[["allen_metadata"]]
ALLEN_MATRIX = snakemake@params[["allen_matrix"]]

OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]
 
STEP2 = "02_seurat/"
STEP3 = "03_sims/"

#-------------------------------------
# Load the metadata of the reference
#-------------------------------------
reference_metadata <- fread(file = file.path(REF, ALLEN_METADATA))
print("Reference metadata loaded. It's size is :")
print(dim(reference_metadata))

########################################## PROPRE AU PROJET ##########################################
#-------------------------------------
# Keep only some part of the cortex
# we want to study: Ssp
#-------------------------------------
reference_metadata <- reference_metadata[reference_metadata$region_label %in% "SSp",]
print("Only SSp are kept in the metadata: The new size of this one is :")
print(dim(reference_metadata))

######################################################################################################

#-------------------------------------
# For the training we need at least
# two cells per labels
#-------------------------------------
occurence <- table(reference_metadata$subclass_label)
one_cell_per_label <- names(occurence[occurence == 1])
reference_metadata <- reference_metadata[!(reference_metadata$subclass_label %in% one_cell_per_label), ]
print("The new size after removing labels thaht only appear once :")
print(dim(reference_metadata))

#-------------------------------------
# Load the matrix of reference
#-------------------------------------
reference_matrix <- h5read(file = file.path(REF, ALLEN_MATRIX), "/data/counts")
reference_genes <- h5read(file = file.path(REF, ALLEN_MATRIX), "/data/gene")
reference_cells <- h5read(file = file.path(REF, ALLEN_MATRIX), "/data/samples")

rownames(reference_matrix) <- as.character(reference_cells)
colnames(reference_matrix) <- as.character(reference_genes)

rm(list = c("reference_genes","reference_cells"))
gc()

print("Reference matrix loaded. The size is :")
print(dim(reference_matrix))

#-------------------------------------
# Keep cells that match both the 
# metadata an matrix of reference
#-------------------------------------
reference_matrix <- reference_matrix[rownames(reference_matrix) %in% reference_metadata$sample_name,]
reference_metadata <- reference_metadata[reference_metadata$sample_name %in% rownames(reference_matrix),]

print("The metadata and matrix of reference have the same cells now : The new size of them tow :")
print(dim(reference_matrix))
print(dim(reference_metadata))

#-------------------------------------
# Order metadata according to matrix
# (Necessary ?)
#-------------------------------------
reference_metadata <- reference_metadata[match(rownames(reference_matrix),reference_metadata$sample_name),]
print("Cells are in the same order in booth files")

#-------------------------------------
# Write the the metadata as csv   
#-------------------------------------
dir.create(file.path(OUTPUTDIR, STEP3, SAMPLE_ID))

fwrite(x = reference_metadata, file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_METADATA, ".csv")))

rm("reference_metadata")
gc()

print("Metadata is written")
#-------------------------------------
# Normalization of the matrix of
# reference
#-------------------------------------
reference_matrix <-  NormalizeData( object = reference_matrix, normalization.method = NORM_METHOD, scale.factor = NORM_SCALE_FACTOR, verbose = FALSE)
print("The matrix of reference is now normalized the same way as our matrix to analyze: Here the first 20 rows and 10 columns :")
print(reference_matrix[1:20, 1:10])

#-------------------------------------
# Load our matrix to annotate
#-------------------------------------
matrix <- fread(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_normalized_matrix.csv"))) # The count matrix is used for this reference.
matrix <- as.matrix(matrix)
rownames(matrix) <- matrix[,"V1"]

print("Our matrix to annotate is loaded, see the dimension below :")
print(dim(matrix))

#-------------------------------------
# Look at genes who are common in both
# matrix and subset genes who are nots
#-------------------------------------
gene_intersection <- intersect(colnames(matrix), colnames(reference_matrix))
reference_matrix <- reference_matrix[, gene_intersection]
matrix <- matrix[, gene_intersection]

print("Only the gene in common in booth our matrix are kept : The new dimmension :")
print(dim(reference_matrix))
print(dim(matrix))

#-------------------------------------
# Wtrite the two new matrix
#-------------------------------------
class(matrix) <- "numeric"
write.csv(matrix, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_MATRIX,".csv")))
print("Our matrix is written")

rm("matrix")
gc()

# reference_matrix <- as.matrix(reference_matrix)
# rownames(reference_matrix) <- reference_matrix[,"sample_name"]
# class(reference_matrix) <- "numeric"
write.csv(reference_matrix, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_MATRIX,".csv")))
print("Our matrix of reference is written")

rm("reference_matrix")
gc()

#-------------------------------------
# create the output file
#-------------------------------------
output_file<-file(TEXT_OUTPUT)
writeLines(c("Files preparation for the Allen reference for SIMS finished (CSV format)"), output_file)
close(output_file)
