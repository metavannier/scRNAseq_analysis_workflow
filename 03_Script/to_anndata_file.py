#..............................................................
# Convert csv matrix to anndata matrix
#..............................................................

#-------------------------------------
# Import
#-------------------------------------
import anndata as an
import os
import sys

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = os.getcwd()
OUTPUTDIR = os.path.join(DIRECTORY, "05_Output")
REF = os.path.join(DIRECTORY, "01_Reference")

SAMPLE_ID = sys.argv[1]
OUTPUT_NAME_REF_MATRIX = sys.argv[2]
OUTPUT_NAME_MATRIX = sys.argv[3]

STEP3 = "03_sims/"
#-------------------------------------
# Write anndata for the reference
# matrix
#-------------------------------------
reference_matrix = an.read_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_REF_MATRIX + ".csv"))
print("reference matrix loaded")
reference_matrix.write_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_REF_MATRIX + ".h5ad"))
print("reference matrix writed")

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
# Write anndata for the matrix
# we want to annotate
#-------------------------------------
matrix = an.read_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_MATRIX + ".csv"))
print("matrix loaded")
matrix.write_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_MATRIX + ".h5ad"))
print("matrix writed")

anndata_matrix = an.read_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, OUTPUT_NAME_MATRIX + ".h5ad"))
print("MATRIX")
print("-------------------------------------------------------------------")
print(anndata_matrix)
print("-------------------------------------------------------------------")
print(anndata_matrix.var_names)
print("-------------------------------------------------------------------")
print(anndata_matrix.obs_names)
print("-------------------------------------------------------------------")