#..............................................................
# Preparation of new metadatas using Remi work for the Allen
# and Di Bella references
#..............................................................

#-------------------------------------
# Library
#-------------------------------------
library(Seurat)
library(data.table)
library(ggplot2)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
REF = file.path((DIRECTORY), "01_Reference")

STEP2 = "02_seurat/"
STEP3 = "03_sims/"
STEP_REMI = "03_remi/"

dir.create(file.path(OUTPUTDIR, STEP_REMI))

#-------------------------------------
# Fichier celltype manipulation
#-------------------------------------

# Open the RDS and write the structure as a .txt file :
remi_celltypes <- readRDS(file.path(REF, "remi/GECortex_PostMsub_celltypefinal.rds" )) # ~ 16 GB
sink(file.path(OUTPUTDIR, STEP_REMI, "remi_celltypefinal_structure.txt" ))
str(remi_celltypes)
sink()

# Write the metadata as a csv file :
remi_celltypes_metadata <- remi_celltypes@meta.data
write.csv(remi_celltypes_metadata, file.path(OUTPUTDIR, STEP_REMI, "remi_celltypefinal_metadata.csv"), row.names = TRUE)

# write metadata for both Arlotta and Allen reference :
motif_arlotta <- "J.Di Bella et al., Nature 2021"
remi_arlotta_metadata <- subset(remi_celltypes_metadata, grepl(motif_arlotta, remi_celltypes_metadata[["study_label"]]))
write.csv(remi_arlotta_metadata, file.path(OUTPUTDIR, STEP_REMI, "remi_arlotta_metadata.csv"), row.names = TRUE)

motif_allen <- "Yao et al., Cell 2021"
remi_allen_metadata <- subset(remi_celltypes_metadata, grepl(motif_allen, remi_celltypes_metadata[["study_label"]]))
remi_allen_metadata <- remi_allen_metadata[remi_allen_metadata$platform_label %in% "10X", ] # I only use this matrix in my analysis
write.csv(remi_allen_metadata, file.path(OUTPUTDIR, STEP_REMI, "remi_allen_metadata.csv"), row.names = TRUE)

#-------------------------------------
# ARLOTTA 
#-------------------------------------
# Open the metadata
arlotta_metadata <- fread(file = file.path(REF, "arlotta/metaData_scDevSC.txt"))

# Remove every cells that are not Non Neuronal
strings_to_remove <- c("Interneurons", "Cajal Retzius cells", "UL CPN", "Layer 4", "DL CPN", "DL_CPN_1", "DL_CPN_2", "SCPN", "NP", "CThPN", "Layer 6b", "Doublet", "group", "Low quality cells", "Red blood cells")
for (pattern in strings_to_remove) {
  rows_to_remove <- grep(pattern, arlotta_metadata$New_cellType)
  arlotta_metadata <- arlotta_metadata[-rows_to_remove, ]
}

# Keep only some column in both files for a futur merge
arlotta_metadata <- subset(arlotta_metadata, select = c("NAME", "New_cellType", "orig_ident")) # Cells, subclass (for NN), and the age of the mouse
remi_arlotta_metadata <- subset(remi_arlotta_metadata, select = c("sample_name", "celltype_label", "biosample_id"))# Cells, celltype (for GABAergic and Glutamatergic), and the age of the mouse

# Give the two data table the same column name for the merge
setnames(arlotta_metadata, "NAME", "sample_name")
setnames(arlotta_metadata, "New_cellType", "celltype_label")
setnames(arlotta_metadata, "orig_ident", "biosample_id")

# Merge the two files
final_arlotta <- rbind(remi_arlotta_metadata, arlotta_metadata)
write.csv(final_arlotta, file.path(OUTPUTDIR, STEP_REMI, "final_arlotta.csv"), row.names = FALSE)

#-------------------------------------
# ALLEN
#-------------------------------------
# Open the metadata
allen_metadata <- fread(file = file.path(REF, "allen/metadata.csv"))

# Keep only cells that are Non Neuronal
strings_to_keep <- c("Astro", "Endo", "Micro-PVM", "Oligo", "SMC-Peri", "VLMC")
allen_metadata <- allen_metadata[allen_metadata$subclass_label %in% strings_to_keep]

# Keep only some column in both files for a futur merge
allen_metadata <- subset(allen_metadata, select = c("sample_name", "subclass_label"))# Cells, subclass (for NN)
remi_allen_metadata <- subset(remi_allen_metadata, select = c("exp_component_name", "celltype_label"))

# Give the two data table the same column name for the merge
setnames(remi_allen_metadata, "exp_component_name", "sample_name")
setnames(allen_metadata, "subclass_label", "celltype_label")

# Merge the two files
final_allen <- rbind(remi_allen_metadata, allen_metadata)
write.csv(final_allen, file.path(OUTPUTDIR, STEP_REMI, "final_allen.csv"), row.names = FALSE)







