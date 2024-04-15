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

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = os.getcwd()
REF = os.path.join(DIRECTORY, "01_Reference")
OUTPUTDIR = os.path.join(DIRECTORY, "05_Output")

REFERENCE_NAME = sys.argv[1]
REFERENCE_MATRIX = sys.argv[2]
PROJECT_NAME = sys.argv[3]
CLASS_LABEL = sys.argv[4]
NUM_WORKERS = sys.argv[5]
MONITOR = sys.argv[6]
PATIENCE = sys.argv[7]
MAX_EPOCH = sys.argv[8]
MATRIX = sys.argv[9]
KEY = sys.argv[10]

STEP3 = "03_sims/"

#-----------------------------------------------------------------------------------------------
#                                       To train the model 
#-----------------------------------------------------------------------------------------------

# -------------------------------------
# Logger : Allows you to visualize 
# data training.
# -------------------------------------
### To use wandb
wandb.login(key = KEY)
logger = WandbLogger(project = PROJECT_NAME, name = PROJECT_NAME, save_dir = os.path.join(OUTPUTDIR, STEP3, REFERENCE_NAME), version = "")

### To turn off wandb
# logger = WandbLogger(offline=True)

#-------------------------------------
# Load the anndata file
#-------------------------------------
reference_matrix = an.read_h5ad(os.path.join(OUTPUTDIR, STEP3, REFERENCE_NAME, REFERENCE_MATRIX + "_join_metadata.h5ad"))

#-------------------------------------
# Custom training jobs
#-------------------------------------
sims = SIMS(data = reference_matrix, class_label = CLASS_LABEL, num_workers = int(NUM_WORKERS), stratify = True) 

# weighting loss inversely proportional by label freq, helps learn rare cell types (recommended)
sims.setup_model(n_a=64, n_d=64, weights=sims.weights)

sims.setup_trainer(
    logger = logger,
    callbacks = [
        EarlyStopping(
            monitor = MONITOR,
            patience = int(PATIENCE),
            verbose = 1,
        ),
        LearningRateMonitor(logging_interval = "epoch"),
        ModelCheckpoint(
            dirpath = os.path.join(OUTPUTDIR, STEP3, REFERENCE_NAME),
            every_n_epochs = 1,
            save_top_k = -1
        ),
    ],
    max_epochs = int(MAX_EPOCH)
)

sims.train()
