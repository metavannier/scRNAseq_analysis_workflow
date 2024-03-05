#..............................................................
# SIMS : Scalable, Interpretable Machine Learning for 
# Single-Cell : A machine learning tool for scRNA-Seq label
# transfer in neuroscience
#..............................................................

#-------------------------------------
# Import
#-------------------------------------
from scsims import SIMS
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning.callbacks import EarlyStopping, LearningRateMonitor, ModelCheckpoint
import os
import sys
import anndata as an
import wandb
# import tensorflow as tf
# print(tf.config.list_physical_devices('GPU'))

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = os.getcwd()
OUTPUTDIR = os.path.join(DIRECTORY, "05_Output")

SAMPLE_ID = sys.argv[1]
MATRIX = sys.argv[2]
EPOCH = sys.argv[3]

STEP3 = "03_sims/"

#-----------------------------------------------------------------------------------------------
#                                To predict labels on our data 
#-----------------------------------------------------------------------------------------------

#-------------------------------------
# Load the model
#-------------------------------------
model = os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, EPOCH)

#-------------------------------------
# Prediction
#-------------------------------------
sims = SIMS(weights_path = model)
cell_prediction = sims.predict(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, MATRIX + ".h5ad"))

#-------------------------------------
# Write prediction in a csv file
#-------------------------------------
matrix = an.read_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, MATRIX + ".h5ad"))

cell_prediction.insert(0,"Cells", matrix.obs_names)

cell_prediction.to_csv(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, MATRIX + "_prediction.csv"), index = False)