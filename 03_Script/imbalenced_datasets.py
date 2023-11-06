#..............................................................
# KNNOR : An oversampling technique for imalanced datasets
#..............................................................

#-------------------------------------
# Import
#-------------------------------------
from knnor import data_augment
from sklearn.preprocessing import LabelEncoder
import pandas as pd
import numpy as np
import sys
import os

#-------------------------------------
# Open metadata and matrix
#-------------------------------------
matrix = pd.read_csv("/scratch/lchabot/scRNAseq_analysis_workflow/05_Output/03_sims/P5_WT/arlotta_reference_matrix_P4_less_labels_knnor.csv")
metadata = pd.read_csv("/scratch/lchabot/scRNAseq_analysis_workflow/05_Output/03_sims/P5_WT/arlotta_reference_metadata_P4_less_labels_knnor.csv")

#-------------------------------------
# Check if cells are in the same order
# in booth file
#-------------------------------------
col_matrix = matrix.iloc[:, 0].tolist()
col_metadata = metadata['NAME'].tolist()
if col_matrix == col_metadata:
    print("Les colonnes sont dans le même ordre.")
else:
    print("Les colonnes ne sont pas dans le même ordre.")  ## Les mettre dans le même ordre

#-------------------------------------
# Add the label column from the 
# metadata at the end of the matrix
#-------------------------------------
col_name = metadata['New_cellType']
matrix['Label'] = col_name

#-------------------------------------
# Put the label into numeric class
#-------------------------------------
label_encoder = LabelEncoder()
matrix['Label'] = label_encoder.fit_transform(matrix['Label'])

# Create an asociate dictionary
label_dictionnary = {index: label for index, label in enumerate(label_encoder.classes_)}
print(label_dictionnary)

# Put the first column as rownames
matrix = matrix.set_index(matrix.columns[0])
matrix.index.name = None

#-------------------------------------
# Convert the matrix to a numpy array
#-------------------------------------
matrix = np.array(matrix)

# X : the matrix, y : the labels (numeric values)
X = matrix[:,:-1]
y = matrix[:,-1]

# -------------------------------------
# Use KNNOR
# -------------------------------------
knnor = data_augment.KNNOR()

# X_new : complete data with augmented datapoints
# y_new : Labels including the augmented ones
# X_aug_min : Just the augmented minority points
# y_aug_min : Labels for only augmented minority points
X_new, y_new, X_aug_min, y_aug_min = knnor.fit_resample(X, y)
y_new = y_new.reshape(-1,1)

# new_data contains the augmented data (Matrix + numeric labels )
new_data = np.append(X_new, y_new, axis=1)


print("")
print("**************************************************************************************************************************************")
print("LABEL AUGMENTE")
print(y_aug_min)
unique_values, counts = np.unique(y_aug_min, return_counts=True)
for value, count in zip(unique_values, counts):
    print(f"{value} est répété {count} fois.")
print("**************************************************************************************************************************************")





# Vous pouvez maintenant utiliser ce dictionnaire pour convertir les valeurs numériques en labels d'origine
# valeur_numerique = 1
# label_d_origine = correspondance_labels[valeur_numerique]
# print(f"La valeur numérique {valeur_numerique} correspond à : {label_d_origine}")

