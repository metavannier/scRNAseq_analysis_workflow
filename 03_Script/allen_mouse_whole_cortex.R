#..............................................................
# Preparation of the folowing reference for sims :
# Paola Arlotta reference : 
#..............................................................

#-------------------------------------
# Library
#-------------------------------------
library(rhdf5)
library(data.table)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]
OUTPUTDIR = file.path(dirname(DIRECTORY), "05_Output")
REF = file.path(dirname(DIRECTORY), "01_Reference")

SAMPLE_ID = snakemake@params[["sample_id"]]

ALLEN_METADATA = snakemake@params[["allen_metadata"]]
ALLEN_MATRIX = snakemake@params[["allen_matrix"]]
ALLEN_GENES = snakemake@params[["allen_genes"]]
ALLEN_CELLS =  snakemake@params[["allen_cells"]]

OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]
 
STEP2 = "02_seurat/"
STEP3 = "03_sims/"

#-------------------------------------
# Load the metadata of the reference
#-------------------------------------
reference_metadata <- fread(file = file.path(REF, ALLEN_METADATA))

########################################## PROPRE AU PROJET ##########################################
#-------------------------------------
# Keep only some part of the cortex
# we want to study: Ssp
#-------------------------------------
reference_metadata <- reference_metadata[reference_metadata$region_label %in% "Ssp",]

######################################################################################################

#-------------------------------------
# Load the matrix of reference
#-------------------------------------
reference_matrix <- h5read(file = file.path(REF, ALLEN_MATRIX), "/data/counts")
reference_genes <- h5read(file = file.path(REF, ALLEN_GENES), "/data/gene")
reference_cells <- h5read(file = file.path(REF, ALLEN_CELLS), "/data/samples")

rownames(reference_matrix) <- as.character(reference_cells)
colnames(reference_matrix) <- as.character(reference_genes)

rm(list = c("reference_genes","reference_cells"))
gc()

#-------------------------------------
# Keep cells thaht match both the 
# metadata an matrix of reference
#-------------------------------------
reference_matrix <- reference_matrix[rownames(reference_matrix) %in% reference_metadata$sample_name,]
reference_metadata <- reference_metadata[reference_metadata$sample_name %in% rownames(reference_matrix),]

#-------------------------------------
# Order metadata according to matrix
# (Necessary ?)
#-------------------------------------
reference_metadata <- reference_metadata[match(rownames(reference_matrix),reference_metadata$sample_name),]

#-------------------------------------
# Wtrite the the metadata as csv   
#-------------------------------------
fwrite(x = reference_metadata, file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_METADATA, ".csv")))

rm("reference_metadata")
gc()

#-------------------------------------
# Load our matrix to annotate
#-------------------------------------
matrix <- fread(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_count_matrix.csv"))) # The count matrix is used for this reference.
matrix <- as.matrix(matrix)
rownames(matrix) <- matrix[,"V1"]

#-------------------------------------
# Look at genes who are common in both
# matrix and subset genes who are nots
#-------------------------------------
gene_intersection <- intersect(colnames(matrix), colnames(reference_matrix))
reference_matrix <- reference_matrix[, gene_intersection]
matrix <- matrix[, gene_intersection]

#-------------------------------------
# Wtrite the two new matrix
#-------------------------------------
dir.create(file.path(OUTPUTDIR, STEP3, SAMPLE_ID))

class(matrix) <- "numeric"
write.csv(matrix, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_MATRIX,".csv")))

rm("matrix")
gc()

class(reference_matrix) <- "numeric"
write.csv(reference_matrix, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_MATRIX,".csv")))

rm("reference_matrix")
gc()

#-------------------------------------
# create the output file
#-------------------------------------
output_file<-file(TEXT_OUTPUT)
writeLines(c("Files preparation for the Arlotta reference for SIMS finished (CSV format)"), output_file)
close(output_file)
