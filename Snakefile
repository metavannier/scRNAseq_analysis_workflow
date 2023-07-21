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

# If using conda environment
container: "docker://condaforge/mambaforge:22.11.1-4"

# ----------------------------------------------
# Load config and sample sheet
# ----------------------------------------------

rawsample = pd.read_table(config["sample"]).set_index(["rawsample"], drop=False)
RAWSAMPLE = expand("{rawsample.rawsample}", rawsample = rawsample.itertuples())
# sample = pd.read_table(config["sample"]).set_index(["name","sample","lane"], drop=False)
# sample_id = pd.read_table(config["sample"]).set_index(["id"], drop=False) 
# sample_name = pd.read_table(config["sample"]).set_index(["name"], drop=False) 

# ----------------------------------------------
# Target rules
# ----------------------------------------------

TYPES = config["run"]["types"].split(',')
LIBRARY = config["run"]["library"].split(',')


SAMPLE = config["fastq"]["sname"].split(',')
NPROJ = config["fastq"]["nproject"]
EXPANSION = config["fastq"]["expansion"]
FEATURES = config["diffexp"]["features"].split(',')
CELLMARKER = config["diffexp"]["cell_marker"].split(',')
WT = config["seurat"]["wt"]
AGRR = config["fastq"]["nproject"]
NUM = config["fastq"]["pair"].split(',')
CLUSTER = config["diffexpsubset"]["cluster"].split(',')

rule all:
	input:
		### fastqc ###
		fastqc_output = expand(OUTPUTDIR + "00_clean/fastqc_output.txt"),
		### multiqc ###
		multiqc_output = expand(OUTPUTDIR + "00_clean/multiqc_output.txt"),
		### ReferenceEnhancer to generate a scRNA-seq optimized transcriptomic reference ###
		### To improve with manual curration of the file overlapping_gene_list
		reference_enhancer_output = expand(OUTPUTDIR + "01_cellranger/reference_enhancer_output.txt"),
		### reference for cellranger ###
		ref_cellranger_output = expand(OUTPUTDIR + "01_cellranger/ref_cellranger_output.txt"),
		### Cell Multiplexing with cellranger multi ###
		multiplexing_output = expand(OUTPUTDIR + "01_cellranger/multiplexing_output.txt"),
		# Cellranger count
		# out_cellranger = expand(OUTPUTDIR + "01_cellranger/{sample}/outs/{sample}_web_summary.html", sample=SAMPLE),
		## If you need to aggregate your data
		# aggrcsv = ROOTDIR + "/aggregation.csv",
		# out_aggregate = expand(OUTPUTDIR + "01_cellranger/{agrr}/outs/aggregate_web_summary.html", agrr=AGRR),
		## If demuxiplexing is used
		# bcf = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.bcf",
		# demuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.best",
		# tabdemuxlet = OUTPUTDIR + "01_cellranger/Mix_MM_lines/outs/demuxlet_Mix_MM_lines.tsv",
		# Seurat
		# count_matrix = OUTPUTDIR + "02_seurat/data_matrix/" + NPROJ + "_count_matrix.csv",
		# data_matrix = OUTPUTDIR + "02_seurat/data_matrix/" + NPROJ + "_data_matrix.csv",
		# scale_data_matrix = OUTPUTDIR + "02_seurat/data_matrix/" + NPROJ + "_scale_data_matrix.csv",
		# seurat_report = OUTPUTDIR + "02_seurat/" + NPROJ + "_seurat_report.html",
		# seurat_object = OUTPUTDIR + "02_seurat/" + NPROJ + "_seurat_object.rds",
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
		# defile_subset_tar = OUTPUTDIR + "04_diffexp_subset/differencial_expression_tests.tar.gz",
		# volcano_subset_tar = OUTPUTDIR + "04_diffexp_subset/volcano_plot.tar.gz",
		# sign_up_subset_tar = OUTPUTDIR + "04_diffexp_subset/up_down_regulated_genes_list.tar.gz",
		# violinplot_subset_tar = OUTPUTDIR + "04_diffexp_subset/violin_plot.tar.gz",
		# clean = OUTPUTDIR + "03_diffexp/clean.txt",
  
# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: ENVDIR + "clean.smk"
include: ENVDIR + "cellranger.smk"
# include: ENVDIR + "demuxlet.smk"
# include: ENVDIR + "seurat.smk"
# include: ENVDIR + "diffexp.smk"
# include: ENVDIR + "diffexp_subset.smk"
# include: ENVDIR + "report.smk"
