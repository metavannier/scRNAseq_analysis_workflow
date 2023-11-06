#..............................................................
# Preparation of the folowing reference for sims :
# Paola Arlotta reference :
# https://www.nature.com/articles/s41586-021-03670-5
#..............................................................

#-------------------------------------
# Library
#-------------------------------------
library(Matrix)
library(data.table)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
REF = file.path((DIRECTORY), "01_Reference")

PATTERN_TO_KEEP = snakemake@params[["pattern_to_keep"]]
SAMPLE_ID = snakemake@params[["sample_id"]]

ARLOTTA_METADATA = snakemake@params[["arlotta_metadata"]]
ARLOTTA_MATRIX = snakemake@params[["arlotta_matrix"]]
ARLOTTA_CELLS = snakemake@params[["arlotta_cells"]]
ARLOTTA_FEATURES = snakemake@params[["arlotta_features"]]

OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]

STEP2 = "02_seurat/"
STEP3 = "03_sims/"

#-------------------------------------
# Load the metadata of the reference
#-------------------------------------
reference_metadata <- fread(file = file.path(REF, ARLOTTA_METADATA))

#-------------------------------------
# Erase some labels in the metadata :
# Low quality cells, Doublet,
# Red blood cells, group
#-------------------------------------
patterns_to_delete <- c("Low quality cells", "Doublet", "Red blood cells", "group") 
rows_to_delete <- apply(sapply(patterns_to_delete, grepl, reference_metadata$New_cellType), 1, any)
reference_metadata <- subset(reference_metadata, !rows_to_delete)

########################################## PROPRE AU PROJET ##########################################
#-------------------------------------
# Rename some subclass to match the
# Allen reference
#-------------------------------------
# reference_metadata$New_cellType[reference_metadata$New_cellType == "UL CPN"] <- "IT"
# reference_metadata$New_cellType[reference_metadata$New_cellType == "Layer 4"] <- "IT"
# reference_metadata$New_cellType[reference_metadata$New_cellType == "DL CPN"] <- "IT"
# reference_metadata$New_cellType[reference_metadata$New_cellType == "DL_CPN_1"] <- "IT"
# reference_metadata$New_cellType[reference_metadata$New_cellType == "DL_CPN_2"] <- "IT"
# reference_metadata$New_cellType[reference_metadata$New_cellType == "SCPN"] <- "L5 PT"
# reference_metadata$New_cellType[reference_metadata$New_cellType == "NP"] <- "L5 NP"
# reference_metadata$New_cellType[reference_metadata$New_cellType == "CThPN"] <- "L6 CT"
# reference_metadata$New_cellType[reference_metadata$New_cellType == "Layer 6b"] <- "L6b"
# reference_metadata$New_cellType[reference_metadata$New_cellType == "Cajal Retzius cells"] <- "CR"

######################################################################################################

#-------------------------------------
# Load the matrix of reference
#-------------------------------------
reference_matrix <- readMM(file = file.path(REF, ARLOTTA_MATRIX))

feature.names = read.delim(file.path(REF, ARLOTTA_FEATURES),
                           header = FALSE,
                           stringsAsFactors = FALSE)

cells.names = read.delim(file.path(REF, ARLOTTA_CELLS),
                           header = FALSE,
                           stringsAsFactors = FALSE)

colnames(reference_matrix) = cells.names$V1
rownames(reference_matrix) = feature.names$V1

# #-------------------------------------
# # Transpose the matrix :
# # Cells in row and genes in columns +
# # Transform rownames and colnames
# # as characters (for sims).
# #-------------------------------------
reference_matrix <- t(reference_matrix)
reference_matrix <- as.matrix(reference_matrix)

colnames(reference_matrix) <- as.character(colnames(reference_matrix)) # Genes
rownames(reference_matrix) <- as.character(rownames(reference_matrix)) # Cells

#-------------------------------------
# Write a new matrix : Keep only
# mouse age that are interesting
# for the project
#-------------------------------------
pattern_to_keep <- c(PATTERN_TO_KEEP) 
rows_to_keep <- grepl(paste(pattern_to_keep, collapse = "|"), rownames(reference_matrix))
reference_matrix <- reference_matrix[rows_to_keep, ]

#-------------------------------------
# Keep cells that match both the 
# metadata an matrix of reference
#-------------------------------------
column_metadata <- reference_metadata$NAME
matching_cells <- intersect(rownames(reference_matrix), column_metadata)
reference_metadata <- reference_metadata[reference_metadata$NAME %in% matching_cells, ]
reference_matrix <- reference_matrix[matching_cells, ]

#-------------------------------------
# Order metadata according to matrix
# (Necessary ?)
#-------------------------------------
reference_metadata <- reference_metadata[match(rownames(reference_matrix),reference_metadata$NAME),]

#-------------------------------------
# Load our matrix to annotate
#-------------------------------------
matrix <- fread(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_normalized_matrix.csv"))) # The normalized matrix is used for this reference.
matrix <- as.matrix(matrix)
rownames(matrix) <- matrix[,"V1"]

#-------------------------------------
# Look at genes who are common in both
# matrix and subset genes who are not
#-------------------------------------
gene_intersection <- intersect(colnames(matrix), colnames(reference_matrix))
reference_matrix <- reference_matrix[, gene_intersection]
matrix <- matrix[, gene_intersection]

#-------------------------------------
# Wtrite the two new matrix +
# create directory
#-------------------------------------
dir.create(file.path(OUTPUTDIR, STEP3, SAMPLE_ID))

class(reference_matrix) <- "numeric"
write.csv(reference_matrix, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_MATRIX,".csv")))

class(matrix) <- "numeric"
write.csv(matrix, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_MATRIX,".csv")))

fwrite(x = reference_metadata, file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_METADATA,".csv")))

#-------------------------------------
# create the output file
#-------------------------------------
output_file<-file(TEXT_OUTPUT)
writeLines(c("Files preparation for the Arlotta reference for SIMS finished (CSV format)"), output_file)
close(output_file)