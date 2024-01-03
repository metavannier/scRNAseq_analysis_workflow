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
# Path to files / folders
#-------------------------------------
DIRECTORY = os.getcwd()
OUTPUTDIR = os.path.join(DIRECTORY, "05_Output")

SAMPLE_ID = sys.argv[1]
MATRIX = sys.argv[2]
METADATA = sys.argv[3]
CLASS_LABEL = sys.argv[4]
MAX_OVERSAMPLING = sys.argv[5]
CELLS_COLUMN = sys.argv[6]

STEP3 = "03_sims/"

#-------------------------------------
# Open metadata and matrix 
#-------------------------------------
print("****************************************************************************************************************************************************************************************")
matrix = pd.read_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, MATRIX + ".csv"))
print("The original matrix :")
print(matrix)
print("")

metadata = pd.read_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, METADATA + ".csv"))
print("The original metadata :")
print(metadata)
print("")

#-------------------------------------
# Occurence of labels in metadata
#-------------------------------------
max_labels = metadata[CLASS_LABEL].value_counts()
print("The data will be increased to the level of the largest occurrence (see table below): ")
print(max_labels)
print("")

#-------------------------------------
# Supress the top occurences in
# metadata if needed !
#-------------------------------------
MAX_OVERSAMPLING = MAX_OVERSAMPLING.split(',')
MAX_OVERSAMPLING = '|'.join(MAX_OVERSAMPLING)

new_metadata = metadata[~metadata[CLASS_LABEL].str.contains(MAX_OVERSAMPLING)]
print("The new metadata after supression of top occurence(s) :")
print(new_metadata)
print("")

#-------------------------------------
# Delete cells in matrix that are not
# anymore in the metadata
#-------------------------------------
new_matrix = matrix[matrix.iloc[:, 0].isin(new_metadata[CELLS_COLUMN])]
print("The new matrix that match the new metadata :")
print(new_matrix)
print("")

#-------------------------------------
# Check if cells are in the same order
# in booth file
#-------------------------------------
are_identical = new_matrix.iloc[:, 0].equals(new_metadata[CELLS_COLUMN])

if are_identical:
    print("Column are in the same order")
else:
    print("Column are not in the same order : it will be reorder") ### A faire !!!


#-------------------------------------
# Add the label column from the 
# metadata at the end of the matrix
#-------------------------------------
new_matrix['Label'] = new_metadata[CLASS_LABEL].tolist()

#-------------------------------------
# Put the label into numeric class
#-------------------------------------
label_encoder = LabelEncoder()
new_matrix['Label'] = label_encoder.fit_transform(new_matrix['Label'])

# Create an asociate dictionary
label_dictionnary = {index: label for index, label in enumerate(label_encoder.classes_)}
print("Labels dictionnary : ")
print(label_dictionnary)
print("")

# Put the first column as rownames
new_matrix = new_matrix.set_index(new_matrix.columns[0])
new_matrix.index.name = None

#-------------------------------------
# Convert the matrix to a numpy array
#-------------------------------------
new_matrix = np.array(new_matrix)

# X : the matrix, y : the labels (numeric values)
X = new_matrix[:,:-1]
y = new_matrix[:,-1]

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

# Transform value numeric back to the labels name
y_aug_min = np.array([label_dictionnary[label] for label in y_aug_min])

# Look at label augmentation
print("Labels augmentation : ")
unique_values, counts = np.unique(y_aug_min, return_counts=True)
for value, count in zip(unique_values, counts):
    print(f"{value} is repeated {count} times.")
print("")

# -------------------------------------
# Create final metadata
# -------------------------------------
# Create a dictionary to store new data
metadata_knnor = {
    CELLS_COLUMN: ["KNNOR" + str(i) for i in range(1, len(y_aug_min) + 1)],
    CLASS_LABEL: list(y_aug_min)
}

# Create a DataFrame from the dictionary
metadata_knnor = pd.DataFrame(metadata_knnor)
print("The knnor metadata only with augmented points :")
print(metadata_knnor)
print("")

# Add the new DataFrame to the metadata
metadata = pd.concat([metadata, metadata_knnor], ignore_index=True)

# Replace all NaN values with -100
metadata.fillna(-100, inplace = True)
print("The metadata with augmented values :")
print(metadata)
print("")

# Occurence of labels with new data created
max_labels = metadata[CLASS_LABEL].value_counts()
print("Occurence of labels with new data created: ")
print(max_labels)
print("")

# Write the new metadata with augmented values
metadata.to_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, "KNNOR_" + METADATA + ".csv"), index = False)

# -------------------------------------
# Create final matrix
# -------------------------------------
# Recover new cell names created for new data
cell_names = list(metadata_knnor[CELLS_COLUMN])

# Add a new column (first column) to the numpy array (for the concat step at the end)
num_rows = X_aug_min.shape[0]
empty_column = np.zeros((num_rows, 1), dtype=int)
X_aug_min = np.hstack((empty_column, X_aug_min))

# Transform the numpy matrix with augmented data to a DataFrame
X_aug_min = pd.DataFrame(X_aug_min, columns = matrix.columns)

# Add the cell names to the first column
X_aug_min.iloc[:, 0] = cell_names
print("The knnor matrix only with augmented points :")
print(X_aug_min)
print("")

# Add the augmented data matrix to our matrix 
matrix = pd.concat([matrix, X_aug_min], ignore_index=True)
print("The matrix with augmented values :")
print(matrix)

# Write the new matrix with augmented data
matrix.to_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, "KNNOR_" + MATRIX + ".csv"), index = False)