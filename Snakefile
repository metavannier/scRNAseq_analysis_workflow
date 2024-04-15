import pandas as pd
import numpy as np
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.1.2")

##### load config and sample sheets #####
configfile: "config.yaml"
validate(config, schema="06_Schemas/config.schema.yaml")

#### Set variables ####
ROOTDIR = os.getcwd()
RAWDATA = srcdir("00_RawData/")
REF = srcdir("01_Reference/")
CONTAINER = srcdir("02_Container/")
SCRIPTDIR = srcdir("03_Script/")
ENVDIR = srcdir("04_Workflow/")
OUTPUTDIR = srcdir("05_Output/")
REPORT = srcdir("07_Report/")
BENCHMARK = srcdir("08_benchmark/")
LOG = srcdir("09_log/")

# If running the workflow with slurm
run_slurm = config["run_slurm"]
if run_slurm:
	container: CONTAINER + "mambaforge:23.1.0-1.sif"

# ----------------------------------------------
# Load sample sheet
# ----------------------------------------------

rawsample = pd.read_table(config["sample"]).set_index(["rawsample"], drop=False)
sample_id = pd.read_table(config["sample"]).set_index(["id"], drop=False) 
sample_name = pd.read_table(config["sample"]).set_index(["name"], drop=False) 

RAWSAMPLE = expand("{rawsample.rawsample}", rawsample = rawsample.itertuples()) # 2_S1_L001, 3_S1_L001 ...
SAMPLE_ID = expand("{sample_id.id}", sample_id = sample_id.itertuples()) # P5_KO,P5_WT,P30_KO,P30_WT
SAMPLE_NAME = expand("{sample_name.name}", sample_name = sample_name.itertuples()) # 2,3,6,7

# ----------------------------------------------
# Target rules
# ----------------------------------------------

TYPES = config["run"]["types"].split(',')
LIBRARY = config["run"]["library"].split(',')

SAMPLE = config["fastq"]["sname"].split(',')
NPROJ = config["fastq"]["nproject"]
FEATURES = config["diffexp"]["features"].split(',')
CELLMARKER = config["diffexp"]["cell_marker"].split(',')
AGRR = config["fastq"]["nproject"] # => A voir car varibale identique à NPROJ
NUM = config["fastq"]["pair"].split(',')
CLUSTER = config["diffexpsubset"]["cluster"].split(',')

rule all:
	input:
		### fastqc ###
		# fastqc_output = expand(OUTPUTDIR + "00_clean/fastqc_output.txt"),
		### multiqc ###
		# multiqc_output = expand(OUTPUTDIR + "00_clean/multiqc_output.txt"),
		# raw_multiqc_html = expand(OUTPUTDIR + "00_clean/raw_multiqc.html"),
		### ReferenceEnhancer to generate a scRNA-seq optimized transcriptomic reference ###
		### To improve with manual curration of the file overlapping_gene_list
		# reference_enhancer_output = expand(OUTPUTDIR + "01_cellranger/reference_enhancer_output.txt"),
		# ### reference for cellranger ###
		# ref_cellranger_output = expand(OUTPUTDIR + "01_cellranger/ref_cellranger_output.txt"),
		### Cell Multiplexing with cellranger multi ###
		# cellranger_output = expand(OUTPUTDIR + "01_cellranger/cellranger_output.txt"),
		# cellranger_html = expand(OUTPUTDIR + "01_cellranger/{sample_id}/{sample_id}_web_summary.html", sample_id = SAMPLE_ID),
		### If you need to aggregate your data
		# aggrcsv = ROOTDIR + "/aggregation.csv",
		# out_aggregate = expand(OUTPUTDIR + "01_cellranger/{agrr}/outs/aggregate_web_summary.html", agrr=AGRR),
		### If demuxiplexing is used
		# bcf = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.bcf",
		# demuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.best",
		# tabdemuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.tsv",
		### Seurat
		# seurat_report = expand(OUTPUTDIR + "02_seurat/{sample_id}/{sample_id}_seurat_report.html", sample_id = SAMPLE_ID),
		### Remi data
		# remi_output = expand(OUTPUTDIR + "remi_output.txt"),
		### Prepare data for SIMS
		# data_for_sims_output = expand(OUTPUTDIR + "03_sims/data_for_sims_output.txt"),
		# anndata_for_sims_output = expand(OUTPUTDIR + "03_sims/anndata_for_sims_output.txt"),
		### KNNOR
		# knnor_output = expand(OUTPUTDIR +"03_sims/knnor_output.txt"),
		# anndata_for_sims_output = expand(OUTPUTDIR + "03_sims/anndata_for_sims_output.txt"),
		### SIMS
		# sims_training_output = expand(OUTPUTDIR + "03_sims/output_sims_training.txt"),
		# sims_prediction_output = expand(OUTPUTDIR + "03_sims/output_sims_prediction.txt"),
		# sims_prediction_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_data_matrix_prediction.csv", sample_id = SAMPLE_ID),
		# sims_prediction_unknown_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_data_matrix_prediction_filtered.csv", sample_id = SAMPLE_ID),
		### Representing cellular assignation on UMAP
		umapAssignation_output = expand(OUTPUTDIR + "03_sims/umapAssignation_output.txt"),
		# umap_sims_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_umap_sims.pdf", sample_id = SAMPLE_ID),
		# umap_sims_threshold_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_umap_sims_threshold.pdf", sample_id = SAMPLE_ID),
		# umap_per_labels_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_umap_per_labels.pdf", sample_id = SAMPLE_ID),

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

run_demultiplex = config["run_demultiplex"]
run_multiplex = config["run_multiplex"]

if run_demultiplex:
	# include: ENVDIR + "clean_demultiplex.smk"
	# include: ENVDIR + "cellranger_demultiplex.smk"
	# include: ENVDIR + "seurat.smk"
	# include: ENVDIR + "remi.smk"
	# include: ENVDIR + "prepare_data_sims.smk"
	# include: ENVDIR + "SIMS.smk"
	# include: ENVDIR + "knnor.smk"
	include: ENVDIR + "umapCellAssignation.smk"

if run_multiplex:	
	# include: ENVDIR + "clean.smk"
	# include: ENVDIR + "cellranger.smk"
	# include: ENVDIR + "seurat.smk"
	# include: ENVDIR + "remi.smk"
	# include: ENVDIR + "prepare_data_sims.smk"
	# include: ENVDIR + "SIMS.smk"
	# include: ENVDIR + "knnor.smk"
	include: ENVDIR + "umapCellAssignation.smk"

