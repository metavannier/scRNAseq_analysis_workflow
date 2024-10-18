#=======================================================================================
# Create the arlotta reference with Remi's work
#=======================================================================================

#-------------------------------------
# Library
#-------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
REF = file.path((DIRECTORY), "01_Reference")

#====================================
# Open the arlotta metadata :
# https://singlecell.broadinstitute.org/single_cell/study/SCP1290/molecular-logic-of-cellular-diversification-in-the-mammalian-cerebral-cortex
#====================================
arlotta_metadata <- fread(file = file.path(REF, "arlotta/metaData_scDevSC.txt"))

### Remove some low quality cells an keep only NN, Migrating/Immature and Dorsal Pallium Progenitor :
patterns_to_delete <- c("Interneurons", "Cajal Retzius cells", "UL CPN", "Layer 4", "DL CPN", "DL_CPN_1", "DL_CPN_2", "SCPN", "NP", "CThPN", "Layer 6b", "Doublet", "group", "Low quality cells", "Red blood cells") 
rows_to_delete <- apply(sapply(patterns_to_delete, grepl, arlotta_metadata$New_cellType), 1, any)
arlotta_metadata <- subset(arlotta_metadata, !rows_to_delete)

### Keep only some column
arlotta_metadata <- subset(arlotta_metadata, select = c("NAME", "New_cellType", "orig_ident")) # Cells, subclass (for NN), and the age of the mouse

### Add the class that corespond to each New Celltype => Rémi's thesis
pattern_class <- c("Apical progenitors" = "Dorsal Pallium Progenitor",
                   "Astrocytes" = "Non Neuronal",
                   "Cycling glial cells" = "Non Neuronal",
                   "Endothelial cells" = "Non Neuronal",
                   "Ependymocytes" = "Non Neuronal",
                   "Immature neurons" = "Immature/Migrating",
                   "Intermediate progenitors" = "Dorsal Pallium Progenitor",
                   "Microglia" = "Non Neuronal",
                   "Migrating neurons" = "Immature/Migrating",
                   "Oligodendrocytes" = "Non Neuronal",
                   "Pericytes" = "Non Neuronal",
                   "VLMC" = "Non Neuronal")

assign_class <- function(cell_type, patterns) {
  match <- patterns[grep(cell_type, names(patterns))]
  if (length(match) > 0) {
    return(match)
  } else {
    return(NA)
  }
}

arlotta_metadata$class_label <- sapply(arlotta_metadata$New_cellType, assign_class, pattern_class)
arlotta_metadata$supertype_label <- arlotta_metadata$New_cellType
arlotta_metadata$celltype_label <- arlotta_metadata$New_cellType

### Rename column for a futur merge
setnames(arlotta_metadata, "NAME", "sample_name")
setnames(arlotta_metadata, "New_cellType", "subclass_label")
setnames(arlotta_metadata, "orig_ident", "biosample_id")

#====================================
# Open Remi's metadata : reanotation
# of the arlotta reference
#====================================
rds_remi <- readRDS(file.path(REF, "remi/GECortex_PostMsub_celltypefinal.rds" ))
rds_metadata_remi <- rds_remi@meta.data

rm("rds_remi")
gc()

### Looking for the arlotta annotation in the whole metadata
motif_arlotta <- "J.Di Bella et al., Nature 2021"
remi_arlotta_metadata <- subset(rds_metadata_remi, grepl(motif_arlotta, rds_metadata_remi[["study_label"]]))

### Only keep some column
remi_arlotta_metadata <- subset(remi_arlotta_metadata, select = c("sample_name", "class_label", "biosample_id", "supertype_label", "subclass_label", "celltype_label"))

### Merge the two metadata 
modify_arlotta_metadata <- rbind(remi_arlotta_metadata, arlotta_metadata)

### Create a new column for cell name, to only keep A,T,C,C (to match the matrix) + write it as csv
modify_arlotta_metadata <- as.data.table(modify_arlotta_metadata)
modify_arlotta_metadata <- modify_arlotta_metadata[, new_sample_name := gsub("[^ATCG]", "", sample_name)]
write.csv(modify_arlotta_metadata, file = file.path(REF, "modify_arlotta_metadata.csv"), row.names = FALSE)

#====================================
# Match our arlotta matrices to 
# the new metadata
# https://singlecell.broadinstitute.org/single_cell/study/SCP1290/molecular-logic-of-cellular-diversification-in-the-mammalian-cerebral-cortex
#====================================

### Path of every h5 matrices in a list
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

list_reference_matrix <- list()

### Read every h5 fie
for (x in names(list_files)){
    filepath <- list_files[[x]]
    list_reference_matrix[[x]] <- Read10X_h5(filepath)
}

### Make sure the cellnames are written in the same way as the metadata
clean_colnames <- function(x){
  current_colnames <- colnames(x)
  modified_colnames <- gsub("[^ATCG]", "", current_colnames)
  colnames(x) <- modified_colnames
  return(x)
}
list_reference_matrix <- lapply(list_reference_matrix, clean_colnames)

print("The diferent matrix for diferent developmental time have been loaded in a list")

### Keep cells that match both the metadata and the matrix
metadata <- data.table()

for (x in names(list_reference_matrix)) {
  print("The developmental time :")
  print(x)
  if (x %in% modify_arlotta_metadata$biosample_id) {
    rows_metadata <- modify_arlotta_metadata[modify_arlotta_metadata$biosample_id == x] # Look in the metadata and keep only the rows that corespond to the matrix (Keep P4 rows for the P4 matrix for exemple)
    cell_metadata <- rows_metadata$new_sample_name # Cells name in the metadata (With just A,T,G,C)
    cell_reference_matrix <- colnames(list_reference_matrix[[x]]) # Cells in the matrix
    matching_cells <- intersect(cell_metadata, cell_reference_matrix) # Look at cells that match in metadata and matrix
    
    # Create a temporary metadata for each matrix with only the cells that match and remove the others
    temporary_metadata <- subset(rows_metadata, new_sample_name %in% matching_cells)
    temporary_metadata$new_sample_name <- paste(x, temporary_metadata$new_sample_name, sep="_") # Add the developmental age in front of cell name

    # Append the new metadata for each matrix to the final one
    metadata <- rbind(metadata, temporary_metadata) 
    
    # Now keep only the matching cells in the matrix and remove the others
    matching_indices <- which(cell_reference_matrix %in% matching_cells)
    list_reference_matrix[[x]] <- list_reference_matrix[[x]][, matching_indices, drop=FALSE]
    colnames(list_reference_matrix[[x]]) <- paste(x, colnames(list_reference_matrix[[x]]), sep = "_") # Add the developmental age in front of cell name
  }
}

print("The final size of the metadata : ")
print(dim(metadata))

reference_matrix <- do.call(cbind, list_reference_matrix)
print("The size of the matrix of reference : ")
print(dim(reference_matrix))

#====================================
# Create a Seurat Object + save it
#====================================
rownames(metadata) <- metadata$new_sample_name
seurat_object <- CreateSeuratObject(counts = reference_matrix, meta.data = metadata)

saveRDS(seurat_object, file = file.path(REF, "modify_arlotta_seurat_object.rds"))


################################################################################################ CREATE FILES TO HAVE ALL CLASS FOR SIMS ANNOTATION

# #=====================================
# # Fichier celltype manipulation
# #=====================================
# # Open the RDS and write the structure as a .txt file :
# remi_celltypes <- readRDS(file.path(REF, "remi/GECortex_PostMsub_celltypefinal.rds" )) # ~ 16 GB
# sink(file.path(OUTPUTDIR, STEP_REMI, "remi_celltypefinal_structure.txt" ))
# str(remi_celltypes)
# sink()

# # Write the metadata as a csv file :
# remi_celltypes_metadata <- remi_celltypes@meta.data
# write.csv(remi_celltypes_metadata, file.path(OUTPUTDIR, STEP_REMI, "remi_celltypefinal_metadata.csv"), row.names = TRUE)

# # write metadata for both Arlotta and Allen reference :
# motif_arlotta <- "J.Di Bella et al., Nature 2021"
# remi_arlotta_metadata <- subset(remi_celltypes_metadata, grepl(motif_arlotta, remi_celltypes_metadata[["study_label"]]))
# write.csv(remi_arlotta_metadata, file.path(OUTPUTDIR, STEP_REMI, "remi_arlotta_metadata.csv"), row.names = TRUE)

# # motif_allen <- "Yao et al., Cell 2021"
# remi_allen_metadata <- subset(remi_celltypes_metadata, grepl(motif_allen, remi_celltypes_metadata[["study_label"]]))
# remi_allen_metadata <- remi_allen_metadata[remi_allen_metadata$platform_label %in% "10X", ] # I only use this matrix in my analysis
# write.csv(remi_allen_metadata, file.path(OUTPUTDIR, STEP_REMI, "remi_allen_metadata.csv"), row.names = TRUE)

# #=====================================
# # ARLOTTA 
# #=====================================
# # Open the metadata
# arlotta_metadata <- fread(file = file.path(REF, "arlotta/metaData_scDevSC.txt"))

# # Remove every cells that are not Non Neuronal
# strings_to_remove <- c("Interneurons", "Cajal Retzius cells", "UL CPN", "Layer 4", "DL CPN", "DL_CPN_1", "DL_CPN_2", "SCPN", "NP", "CThPN", "Layer 6b", "Doublet", "group", "Low quality cells", "Red blood cells")
# for (pattern in strings_to_remove) {
#   rows_to_remove <- grep(pattern, arlotta_metadata$New_cellType)
#   arlotta_metadata <- arlotta_metadata[-rows_to_remove, ]
# }

# # Keep only some column in both files for a futur merge
# arlotta_metadata <- subset(arlotta_metadata, select = c("NAME", "New_cellType", "orig_ident")) # Cells, subclass (for NN), and the age of the mouse
# remi_arlotta_metadata <- subset(remi_arlotta_metadata, select = c("sample_name", "supertype_label", "biosample_id"))# Cells, celltype (for GABAergic and Glutamatergic), and the age of the mouse

# # Give the two data table the same column name for the merge
# setnames(arlotta_metadata, "NAME", "sample_name")
# setnames(arlotta_metadata, "New_cellType", "supertype_label")
# setnames(arlotta_metadata, "orig_ident", "biosample_id")

# # Merge the two files
# final_arlotta <- rbind(remi_arlotta_metadata, arlotta_metadata)
# write.csv(final_arlotta, file.path(OUTPUTDIR, STEP_REMI, "supertype_final_arlotta.csv"), row.names = FALSE)

# #=====================================
# # ALLEN
# #=====================================
# # Open the metadata
# allen_metadata <- fread(file = file.path(REF, "allen/metadata.csv"))

# # Keep only cells that are Non Neuronal
# strings_to_keep <- c("Astro", "Endo", "Micro-PVM", "Oligo", "SMC-Peri", "VLMC")
# allen_metadata <- allen_metadata[allen_metadata$subclass_label %in% strings_to_keep]

# # Keep only some column in both files for a futur merge
# allen_metadata <- subset(allen_metadata, select = c("sample_name", "subclass_label"))# Cells, subclass (for NN)
# remi_allen_metadata <- subset(remi_allen_metadata, select = c("exp_component_name", "celltype_label"))

# # Give the two data table the same column name for the merge
# setnames(remi_allen_metadata, "exp_component_name", "sample_name")
# setnames(allen_metadata, "subclass_label", "celltype_label")

# # Merge the two files
# final_allen <- rbind(remi_allen_metadata, allen_metadata)
# write.csv(final_allen, file.path(OUTPUTDIR, STEP_REMI, "final_allen.csv"), row.names = FALSE)

################################################################################################ CREATE FILES TO HAVE REMI'S ANNOTATION 

# remi <- readRDS(file.path(REF, "remi/GECortex_PostMsub_celltypefinal.rds" )) # ~ 16 GB
# cell_names <- grep("^P30_wt_", colnames(remi), value = TRUE)
# remi <- subset(remi, cells = cell_names) 
# saveRDS(remi, file = file.path(OUTPUTDIR, STEP_REMI, "P30_WT_remi.rds"))

################################################################################################ LABEL TRANSFERT FROM THE REFERENCE : ARLOTTA AND ALLEN 