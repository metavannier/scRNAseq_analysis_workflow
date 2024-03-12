#..............................................................
# Preparation of the folowing reference for sims :
# Paola Arlotta reference :
# https://www.nature.com/articles/s41586-021-03670-5
#..............................................................
###########################################################################
#-------------------------------------
# Library
#-------------------------------------
library(Matrix)
library(data.table)
library(Seurat)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
REF = file.path((DIRECTORY), "01_Reference")

DEVELOPMENTAL_TIME = snakemake@params[["developmental_time"]]
SAMPLE_ID = snakemake@params[["sample_id"]]

NORM_METHOD = snakemake@params[["norm_method"]]
NORM_SCALE_FACTOR = snakemake@params[["norm_scale_factor"]]

OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]

STEP2 = "02_seurat/"
STEP3 = "03_sims/"

#-------------------------------------
# Load the metadata of the reference
#-------------------------------------
reference_metadata <- fread(file = file.path(REF, "arlotta/metaData_scDevSC.txt"))
reference_metadata <- reference_metadata[-1, ]
reference_metadata <- reference_metadata[, NEW_NAME := gsub("[^ATCG]", "", NAME)] ### Keep only the A,T,G and C in the original name

#-------------------------------------
# For the training we need at least
# two cells per labels
#-------------------------------------
occurence <- table(reference_metadata$New_cellType)
one_cell_per_label <- names(occurence[occurence == 1])
reference_metadata <- reference_metadata[!(reference_metadata$New_cellType %in% one_cell_per_label), ]

#-------------------------------------
# Erase some labels in the metadata :
# Low quality cells, Doublet,
# Red blood cells, group
#-------------------------------------
patterns_to_delete <- c("Low quality cells", "Doublet", "Red blood cells", "group") 
rows_to_delete <- apply(sapply(patterns_to_delete, grepl, reference_metadata$New_cellType), 1, any)
reference_metadata <- subset(reference_metadata, !rows_to_delete)

########################################################################################################################

### Add the class_label
reference_metadata$class_label <- NA

# Dorsal Pallium Progenitor
reference_metadata$class_label <- ifelse(reference_metadata$New_cellType %in% c("Apical progenitors", "Intermediate progenitors"), "Dorsal Pallium Progenitor", reference_metadata$class_label)

# GABAergic
reference_metadata$class_label <- ifelse(reference_metadata$New_cellType %in% c("Interneurons"), "GABAergic", reference_metadata$class_label)

# Immature_Migrating
reference_metadata$class_label <- ifelse(reference_metadata$New_cellType %in% c("Immature neurons", "Migrating neurons"), "Immature/Migrating", reference_metadata$class_label)

# Glutamatergic
reference_metadata$class_label <- ifelse(reference_metadata$New_cellType %in% c("Cajal Retzius cells", "CThPN", "DL CPN", "DL_CPN_2", "Layer 4", "Layer 6b", "NP", "SCPN", "UL CPN"), "Glutamatergic", reference_metadata$class_label)

# Non Neuronal
reference_metadata$class_label <- ifelse(reference_metadata$New_cellType %in% c("Astrocytes", "Cycling glial cells", "Endothelial cells", "Microglia", "Oligodendrocytes", "Pericytes", "VLMC", "Ependymocytes"), "Non Neuronal", reference_metadata$class_label)

print(table(reference_metadata$class_label))

########################################################################################################################

#-------------------------------------
# Load the matrix of reference 
# (Lot of developmental time points)
#-------------------------------------
# Create a list to store al the file path
list_files <- list(
E10 <- file.path(REF, "arlotta/GSM5277843_E10_v1_filtered_feature_bc_matrix.h5"),
E11 <- file.path(REF, "arlotta/GSM4635072_E11_5_filtered_gene_bc_matrices_h5.h5"),
E12 <- file.path(REF, "arlotta/GSM4635073_E12_5_filtered_gene_bc_matrices_h5.h5"),
E13 <- file.path(REF, "arlotta/GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5"),
E14 <- file.path(REF, "arlotta/GSM4635075_E14_5_filtered_gene_bc_matrices_h5.h5"),
E15 <- file.path(REF, "arlotta/GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5"),
E16 <- file.path(REF, "arlotta/GSM4635077_E16_filtered_gene_bc_matrices_h5.h5"),
E17 <- file.path(REF, "arlotta/GSM5277844_E17_5_filtered_feature_bc_matrix.h5"),
E18_S1 <- file.path(REF, "arlotta/GSM4635078_E18_5_S1_filtered_gene_bc_matrices_h5.h5"),
E18_S3 <- file.path(REF, "arlotta/GSM4635079_E18_S3_filtered_gene_bc_matrices_h5.h5"),
P1 <- file.path(REF, "arlotta/GSM4635081_P1_S2_filtered_gene_bc_matrices_h5.h5"),
P1_S1 <- file.path(REF, "arlotta/GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5"),
P4 <- file.path(REF, "arlotta/GSM5277845_P4_filtered_feature_bc_matrix.h5")
)

list_files <- setNames(list_files, c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18_S1", "E18_S3", "P1", "P1_S1", "P4"))

# file we want to open (depending on what is written in the config.yaml)
file_to_open <- DEVELOPMENTAL_TIME

list_reference_matrix <- list()

for (x in file_to_open){
    filepath <- list_files[names(list_files) == x]
    list_reference_matrix[[x]] <- Read10X_h5(filepath[[x]])
}

#-------------------------------------
# Normalization + Transpose them
# (cell in rows) + Transform them
# in matrix
#-------------------------------------
# list_reference_matrix <- lapply(list_reference_matrix, function(x){
#     normalized_matrix <- NormalizeData( object = x, normalization.method = NORM_METHOD, scale.factor = NORM_SCALE_FACTOR, verbose = FALSE)
# })

list_reference_matrix <-  lapply(list_reference_matrix, function(x) as.matrix(x))
list_reference_matrix <-  lapply(list_reference_matrix, function(x) t(x))

#-------------------------------------
# Keep only A,T,G and C in every 
# matrix for the cells name (to match
# the metadata)
#-------------------------------------
clean_rownames <- function(x){
  current_rownames <- rownames(x)
  modified_rownames <- gsub("[^ATCG]", "", current_rownames)
  rownames(x) <- modified_rownames
  return(x)
}
list_reference_matrix <- lapply(list_reference_matrix, clean_rownames)

#-------------------------------------
# Keep cells that match both the 
# metadata an matrix of reference
#-------------------------------------
metadata <- data.table()

for(x in names(list_reference_matrix)){
    if(x %in% reference_metadata$orig_ident){
        rows_metadata <- reference_metadata[orig_ident == x] # Look in the metadata and keep only the rows that corespond to the matrix (Keep P4 rows for the P4 matrix for exemple)
        cell_metadata <- rows_metadata$NEW_NAME # Cells in the metadata
        cell_reference_matrix <- rownames(list_reference_matrix[[x]]) # Cells in the metadata
        matching_cells <- intersect(cell_metadata, cell_reference_matrix) # Look at cells that match
        temporary_metadata <- subset(rows_metadata, NEW_NAME %in% matching_cells) # Create a temporary metadata for each matrix whith only the cells that match and remove the others
        temporary_metadata$NEW_NAME <- paste(x, temporary_metadata$NEW_NAME, sep="_") # Add the develomental age in front of cell name (some diferent developmental point have same name for cells, that is why)
        metadata <- rbind(metadata, temporary_metadata) # Append the for each matrix the new metadata to the final one
        list_reference_matrix[[x]] <- list_reference_matrix[[x]][matching_cells, ] # Now keep only the mathing cell in the matrix and remove the others
        rownames(list_reference_matrix[[x]]) <- paste(x, rownames(list_reference_matrix[[x]]), sep = "_") # Add the develomental age in front of cell name (Like this also names still match metadata)
    }
}

reference_matrix <- do.call(rbind, list_reference_matrix)

#-------------------------------------
# Load our matrix to annotate
#-------------------------------------
matrix <- fread(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_count_matrix.csv"))) # The normalized matrix is used for this reference.
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

fwrite(x = metadata, file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_METADATA,".csv")))

#-------------------------------------
# create the output file
#-------------------------------------
output_file<-file(TEXT_OUTPUT)
writeLines(c("Files preparation for the Arlotta reference for SIMS finished (CSV format)"), output_file)
close(output_file)

######################################################################################################################################################################################################################################################
# Old code

# #-------------------------------------
# # Path to files / folders
# #-------------------------------------
# DIRECTORY = getwd()
# TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]
# OUTPUTDIR = file.path((DIRECTORY), "05_Output")
# REF = file.path((DIRECTORY), "01_Reference")

# PATTERN_TO_KEEP = snakemake@params[["pattern_to_keep"]]
# SAMPLE_ID = snakemake@params[["sample_id"]]

# ARLOTTA_METADATA = snakemake@params[["arlotta_metadata"]]
# ARLOTTA_MATRIX = snakemake@params[["arlotta_matrix"]]
# ARLOTTA_CELLS = snakemake@params[["arlotta_cells"]]
# ARLOTTA_FEATURES = snakemake@params[["arlotta_features"]]

# OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
# OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
# OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]

# STEP2 = "02_seurat/"
# STEP3 = "03_sims/"

# #-------------------------------------
# # Load the metadata of the reference
# #-------------------------------------
# reference_metadata <- fread(file = file.path(REF, ARLOTTA_METADATA))

# #-------------------------------------
# # Erase some labels in the metadata :
# # Low quality cells, Doublet,
# # Red blood cells, group
# #-------------------------------------
# patterns_to_delete <- c("Low quality cells", "Doublet", "Red blood cells", "group") 
# rows_to_delete <- apply(sapply(patterns_to_delete, grepl, reference_metadata$New_cellType), 1, any)
# reference_metadata <- subset(reference_metadata, !rows_to_delete)

# #-------------------------------------
# # Load the matrix of reference
# #-------------------------------------
# reference_matrix <- readMM(file = file.path(REF, ARLOTTA_MATRIX))

# feature.names = read.delim(file.path(REF, ARLOTTA_FEATURES),
#                            header = FALSE,
#                            stringsAsFactors = FALSE)

# cells.names = read.delim(file.path(REF, ARLOTTA_CELLS),
#                            header = FALSE,
#                            stringsAsFactors = FALSE)

# colnames(reference_matrix) = cells.names$V1
# rownames(reference_matrix) = feature.names$V1

# # #-------------------------------------
# # # Transpose the matrix :
# # # Cells in row and genes in columns +
# # # Transform rownames and colnames
# # # as characters (for sims).
# # #-------------------------------------
# reference_matrix <- t(reference_matrix)
# reference_matrix <- as.matrix(reference_matrix)

# colnames(reference_matrix) <- as.character(colnames(reference_matrix)) # Genes
# rownames(reference_matrix) <- as.character(rownames(reference_matrix)) # Cells

# #-------------------------------------
# # Write a new matrix : Keep only
# # mouse age that are interesting
# # for the project
# #-------------------------------------
# pattern_to_keep <- c(PATTERN_TO_KEEP) 
# rows_to_keep <- grepl(paste(pattern_to_keep, collapse = "|"), rownames(reference_matrix))
# reference_matrix <- reference_matrix[rows_to_keep, ]

# #-------------------------------------
# # Keep cells that match both the 
# # metadata an matrix of reference
# #-------------------------------------
# column_metadata <- reference_metadata$NAME
# matching_cells <- intersect(rownames(reference_matrix), column_metadata)
# reference_metadata <- reference_metadata[reference_metadata$NAME %in% matching_cells, ]
# reference_matrix <- reference_matrix[matching_cells, ]

# #-------------------------------------
# # Order metadata according to matrix
# # (Necessary ?)
# #-------------------------------------
# reference_metadata <- reference_metadata[match(rownames(reference_matrix),reference_metadata$NAME),]

# #-------------------------------------
# # Load our matrix to annotate
# #-------------------------------------
# matrix <- fread(file = file.path(OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_normalized_matrix.csv"))) # The normalized matrix is used for this reference.
# matrix <- as.matrix(matrix)
# rownames(matrix) <- matrix[,"V1"]

# #-------------------------------------
# # Look at genes who are common in both
# # matrix and subset genes who are not
# #-------------------------------------
# gene_intersection <- intersect(colnames(matrix), colnames(reference_matrix))
# reference_matrix <- reference_matrix[, gene_intersection]
# matrix <- matrix[, gene_intersection]

# #-------------------------------------
# # Wtrite the two new matrix +
# # create directory
# #-------------------------------------
# dir.create(file.path(OUTPUTDIR, STEP3, SAMPLE_ID))

# class(reference_matrix) <- "numeric"
# write.csv(reference_matrix, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_MATRIX,".csv")))

# class(matrix) <- "numeric"
# write.csv(matrix, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_MATRIX,".csv")))

# fwrite(x = reference_metadata, file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_METADATA,".csv")))

# #-------------------------------------
# # create the output file
# #-------------------------------------
# output_file<-file(TEXT_OUTPUT)
# writeLines(c("Files preparation for the Arlotta reference for SIMS finished (CSV format)"), output_file)
# close(output_file)