#..............................................................
# 
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

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = os.getcwd()
REF = os.path.join(DIRECTORY, "01_Reference")
OUTPUTDIR = os.path.join(DIRECTORY, "05_Output")

SAMPLE_ID = sys.argv[1]
REFERENCE_MATRIX = sys.argv[2]
PROJECT_NAME = sys.argv[3]
CLASS_LABEL = sys.argv[4]
NUM_WORKERS = sys.argv[5]
MONITOR = sys.argv[6]
PATIENCE = sys.argv[7]
MAX_EPOCH = sys.argv[8]
MATRIX = sys.argv[9]

STEP3 = "03_sims/"

# #-----------------------------------------------------------------------------------------------
# #                                       To train the model 
# #-----------------------------------------------------------------------------------------------

# #-------------------------------------
# # Logger : Allows you to visualize 
# # data training.
# #-------------------------------------
# logger = WandbLogger(project = PROJECT_NAME, name = PROJECT_NAME, save_dir = os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID), version = "")

# #-------------------------------------
# # Load the anndata file
# #-------------------------------------
# reference_matrix = an.read_h5ad(os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, REFERENCE_MATRIX + "_join_metadata.h5ad"))

# #-------------------------------------
# # Custom training jobs
# #-------------------------------------
# sims = SIMS(data = reference_matrix, class_label = CLASS_LABEL, num_workers = 12, stratify = True) 

# # weighting loss inversely proportional by label freq, helps learn rare cell types (recommended)
# sims.setup_model(n_a=64, n_d=64, weights=sims.weights)

# sims.setup_trainer(
#     logger = logger,
#     callbacks = [
#         EarlyStopping(
#             monitor = MONITOR,
#             patience = 25,
#             verbose = 1
#         ),
#         LearningRateMonitor(logging_interval = "epoch"),
#         ModelCheckpoint(
#             dirpath = os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID),
#             filename = "checkpoint"
#         ),
#     ],
#     max_epochs = 100,
#     devices = 1
# )

# sims.train()

#-----------------------------------------------------------------------------------------------
#                                To predict labels on our data 
#-----------------------------------------------------------------------------------------------

#-------------------------------------
# Load the model
#-------------------------------------
model = os.path.join(OUTPUTDIR, STEP3, SAMPLE_ID, "checkpoint-v1.ckpt")

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