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
OUTPUTDIR = file.path(dirname(DIRECTORY), "05_Output")
REF = file.path(dirname(DIRECTORY), "01_Reference")

STEP2 = "02_seurat/"
STEP3 = "03_sims/"

#-------------------------------------
# Load the metadata of the reference                              ### A FINIR
#-------------------------------------
reference_metadata <- fread(file = file.path(REF, "allen_mouse_whole_cortex", ""))

#-------------------------------------
# Load the matrix of reference
#-------------------------------------
reference_matrix <- h5read("")
reference_genes <- h5read("")
reference_cells <- h5read("")

rownames(reference_matrix) <- as.character(reference_cells)
colnames(reference_matrix) <- as.character(reference_genes)

rm(list = c("reference_genes","reference_cells"))
gc()

#-------------------------------------
# Keep cells thaht match both the 
# metadata an matrix of reference
#-------------------------------------


#-------------------------------------
# Order metadata according to matrix
# (Necessary ?)
#-------------------------------------


#-------------------------------------
# Load our matrix to annotate
#-------------------------------------


#-------------------------------------
# Look at genes who are common in both
# matrix and subset genes who are nots
#-------------------------------------


#-------------------------------------
# Wtrite the two new matrix and the
# metadata as csv   
#-------------------------------------