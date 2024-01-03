#..............................................................
# Convert csv matrix to anndata matrix
#..............................................................

#-------------------------------------
# Import
#-------------------------------------
import anndata as an
import pandas as pd
import os
import sys

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = os.getcwd()
REF = os.path.join(DIRECTORY, "01_Reference")
OUTPUTDIR = os.path.join(DIRECTORY, "05_Output")

SAMPLE_ID = sys.argv[1]
OUTPUT_NAME_REF_MATRIX = sys.argv[2]
OUTPUT_NAME_MATRIX = sys.argv[3]
OUTPUT_NAME_REF_METADATA = sys.argv[4]
CELLS_COLUMN = sys.argv[5]

STEP3 = "03_sims/"

#-------------------------------------
# Load reference matrix
#-------------------------------------

reference_matrix = an.read_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_REF_MATRIX + ".csv"))

#-------------------------------------
# Write matrix as anndata
#-------------------------------------
reference_matrix.write_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_REF_MATRIX + ".h5ad"))

anndata_ref_matrix = an.read_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_REF_MATRIX + ".h5ad"))
print("REFERENCE MATRIX")
print("-------------------------------------------------------------------")
print(anndata_ref_matrix)
print("-------------------------------------------------------------------")
print(anndata_ref_matrix.var_names)
print("-------------------------------------------------------------------")
print(anndata_ref_matrix.obs_names)
print("-------------------------------------------------------------------")

#-------------------------------------
# Write matrix as anndata with
# metadata inside
#-------------------------------------
reference_metadata = pd.read_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_REF_METADATA + ".csv"), index_col = CELLS_COLUMN )
anndata_ref_matrix.obs = pd.concat([anndata_ref_matrix.obs, reference_metadata], axis=1, join='inner')

anndata_ref_matrix.write_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_REF_MATRIX + "_join_metadata.h5ad"))

print("REFERENCE MATRIX AVEC METADATA")
print("-------------------------------------------------------------------")
print(anndata_ref_matrix)
print("-------------------------------------------------------------------")
print(anndata_ref_matrix.var_names)
print("-------------------------------------------------------------------")
print(anndata_ref_matrix.obs_names)
print("-------------------------------------------------------------------")

#-------------------------------------
# Load the matrix to annotate
#-------------------------------------
matrix = an.read_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_MATRIX + ".csv"))

#-------------------------------------
# Write anndata for the matrix
# we want to annotate
#-------------------------------------
matrix.write_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_MATRIX + ".h5ad"))

anndata_matrix = an.read_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_MATRIX + ".h5ad"))
print("")
print("")
print("MATRIX")
print("-------------------------------------------------------------------")
print(anndata_matrix)
print("-------------------------------------------------------------------")
print(anndata_matrix.var_names)
print("-------------------------------------------------------------------")
print(anndata_matrix.obs_names)
print("-------------------------------------------------------------------")
