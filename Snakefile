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
		### reference for cellranger ###
		# ref_cellranger_output = expand(OUTPUTDIR + "01_cellranger/ref_cellranger_output.txt"),
		### Cell Multiplexing with cellranger multi ###
		# cellranger_output = expand(OUTPUTDIR + "01_cellranger/cellranger_output.txt"),
		# cellranger_html = expand(OUTPUTDIR + "01_cellranger/{sample_id}/web_summary.html", sample_id = SAMPLE_ID),
		## If you need to aggregate your data
		# aggrcsv = ROOTDIR + "/aggregation.csv",
		# out_aggregate = expand(OUTPUTDIR + "01_cellranger/{agrr}/outs/aggregate_web_summary.html", agrr=AGRR),
		## If demultiplexing is used
		# bcf = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.bcf",
		# demuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.best",
		# tabdemuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.tsv",
		### Seurat ###
		# seurat_report = expand(OUTPUTDIR + "02_seurat/{sample_id}/{sample_id}_seurat_report.html", sample_id = SAMPLE_ID),
		### Normalisation with other method than Seurat ###
        # normalisation_output = expand(OUTPUTDIR + "02_seurat/normalisation_output.txt"),
		# NE MARCHE PAS : Doublet detection
		# doublets_output = expand(OUTPUTDIR + "rm_doublet/doublets_output.txt"),
		### Prepare data for SIMS
		# data_for_sims_output = expand(OUTPUTDIR + "03_sims/data_for_sims_output.txt"),
		# anndata_for_sims_output = expand(OUTPUTDIR + "03_sims/anndata_for_sims_output.txt"),
		### KNNOR
		# knnor_output = expand(OUTPUTDIR +"03_sims/knnor_output.txt"),
		### SIMS
		# sims_training_output = expand(OUTPUTDIR + "03_sims/output_sims_training.txt"),
		# sims_prediction_output = expand(OUTPUTDIR + "03_sims/output_sims_prediction.txt"),
		# sims_prediction_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_data_matrix_prediction.csv", sample_id = SAMPLE_ID),
		# sims_prediction_unknown_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_data_matrix_prediction_filtered.csv", sample_id = SAMPLE_ID),
		unknown_prediction_output = expand(OUTPUTDIR + "03_sims/unknown_prediction_output.txt"),
		evaluate_prediction_output = expand(OUTPUTDIR + "03_sims/evaluate_prediction_output.txt"),
		### Representing cellular assignation on UMAP
		# umapAssignation_output = expand(OUTPUTDIR + "03_sims/umapAssignation_output.txt"),
		# umap_sims_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_umap_sims.pdf", sample_id = SAMPLE_ID),
		# umap_sims_threshold_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_umap_sims_threshold.pdf", sample_id = SAMPLE_ID),
		# umap_per_labels_report = expand(OUTPUTDIR + "03_sims/{sample_id}/{sample_id}_umap_per_labels.pdf", sample_id = SAMPLE_ID),
		## Differential expression analyses
		# violinplot = expand(OUTPUTDIR + "03_diffexp/violin_plot/{features}_violin_plot.pdf", features=FEATURES),
		# umapfeature = expand(OUTPUTDIR + "03_diffexp/umap_plot/{features}_umapfeature_plot.pdf", features=FEATURES),
		# tsnefeature = expand(OUTPUTDIR + "03_diffexp/tsne_plot/{features}_tsnefeature_plot.pdf", features=FEATURES),
		# ridgefeature = expand(OUTPUTDIR + "03_diffexp/ridge_plot/{features}_ridgefeature_plot.pdf", features=FEATURES),
		# heatmapfeature = OUTPUTDIR + "03_diffexp/heatmap/heatmapfeature.pdf",
		# diffexp_report = OUTPUTDIR + "03_diffexp/" + NPROJ + "_diffexp_report.html",
		# defile_allcells = expand(OUTPUTDIR + "03_diffexp/differential_expression_features/{wt}_vs_AllCellMarker_DE.csv", wt=WT),
		# defile = expand(OUTPUTDIR + "03_diffexp/differential_expression_features/{wt}_vs_{cellmarker}_DE.csv", wt=WT, cellmarker=CELLMARKER),
		# volcano_allcluster = expand(OUTPUTDIR + "03_diffexp/volcanoplot/volcano_{wt}_vs_AllCellMarker.pdf", wt=WT),
		# volcano = expand(OUTPUTDIR + "03_diffexp/volcanoplot/volcano_{wt}_vs_{cellmarker}.pdf", wt=WT, cellmarker=CELLMARKER),
		# sign_up_allcell = expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_AllCellMarker_signif-up-regulated.txt", wt=WT),
		# sign_up = expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_{cellmarker}_signif-up-regulated.txt", wt=WT, cellmarker=CELLMARKER),
		# sign_down_allcell = expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_AllCellMarker_signif-down-regulated.txt", wt=WT),
		# sign_down = expand(OUTPUTDIR + "03_diffexp/de_significant/{wt}_vs_{cellmarker}_signif-down-regulated.txt", wt=WT, cellmarker=CELLMARKER),
		### Differential expression analyses inside each cluster after subset (by gene expression)
		# diffexp_subset_report = OUTPUTDIR + "04_diffexp_subset/" + NPROJ + "_diffexp_subset_report.html",
		# defile_subset = expand(OUTPUTDIR + "04_diffexp_subset/differential_expression/{cluster}_subcluster_DE.csv", cluster=CLUSTER),
		# volcano_subset = expand(OUTPUTDIR + "04_diffexp_subset/volcanoplot/volcano_{cluster}.pdf", cluster=CLUSTER),
		# sign_up_subset = expand(OUTPUTDIR + "04_diffexp_subset/de_significant/{cluster}_signif-up-regulated.txt", cluster=CLUSTER),
		# sign_down_subset = expand(OUTPUTDIR + "04_diffexp_subset/de_significant/{cluster}_signif-down-regulated.txt", cluster=CLUSTER),
		# violinplot_subset = expand(OUTPUTDIR + "04_diffexp_subset/violin_plot/{features}_violin_plot.pdf", features=FEATURES),
		### Compress result files and figures ###
		# data_matrix_tar = OUTPUTDIR + "02_seurat/data_matrix.tar.gz",
		# violinplot_tar = OUTPUTDIR + "03_diffexp/violin_plot.tar.gz",
		# umapfeature_tar = OUTPUTDIR + "03_diffexp/umapfeature_plot.tar.gz",
		# tsnefeature_tar = OUTPUTDIR + "03_diffexp/tsnefeature_plot.tar.gz",
		# ridgefeature_tar = OUTPUTDIR + "03_diffexp/ridgefeature_plot.tar.gz",
		# heatmapfeature_tar = OUTPUTDIR + "03_diffexp/heatmapfeature.tar.gz",
		# defile_tar = OUTPUTDIR + "03_diffexp/differencial_expression_tests.tar.gz",
		# volcano_tar = OUTPUTDIR + "03_diffexp/volcano_plot.tar.gz",
		# sign_up_down_tar = OUTPUTDIR + "03_diffexp/up_down_regulated_genes_list.tar.gz",
		# defile_subset_tar = OUTPUTDIR + "04_diffexp_sub		

# ----------------------------------------------
# setup report
# ----------------------------------------------

report: "07_Report/workflow.rst"

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

run_demultiplex = config["run_demultiplex"]
run_multiplex = config["run_multiplex"]

if run_demultiplex:
	include: ENVDIR + "clean_demultiplex.smk"
	include: ENVDIR + "cellranger_demultiplex.smk"
	include: ENVDIR + "seurat.smk"
	include: ENVDIR + "prepare_data_sims.smk"
	include: ENVDIR + "SIMS.smk"
	include: ENVDIR + "umapCellAssignation.smk"
	include: ENVDIR + "knnor.smk"

if run_multiplex:	
	# include: ENVDIR + "clean.smk"
	# include: ENVDIR + "cellranger.smk"
	# NE MARCHE PAS include: ENVDIR + "processing_multiplex.smk"
	# include: ENVDIR + "seurat.smk"
	# NE MARCHE PAS include: ENVDIR + "rm_doublets.smk"
	# include: ENVDIR + "prepare_data_sims.smk"
	include: ENVDIR + "SIMS.smk"
	# include: ENVDIR + "umapCellAssignation.smk"

# include: ENVDIR + "demuxlet.smk"
# include: ENVDIR + "seurat.smk"
# include: ENVDIR + "diffexp.smk"
# include: ENVDIR + "diffexp_subset.smk"
# include: ENVDIR + "report.smk"
