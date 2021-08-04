import pandas as pd
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

# Use glob statement to find all samples in 'raw_data' directory
## Wildcard '{num}' must be equivalent to 'R1' or '1', meaning the read pair designation.
SAMPLE_LIST,NUMS = glob_wildcards(RAWDATA + "/{sample}_L001_{num}_001.fastq")
# # Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
SET_NUMS = set(NUMS)

rule all:
  input:
    # # fastqc
    # raw_html = expand("{outputdir}00_fastqc/{sample}_L001_{num}_001_fastqc.html", outputdir=OUTPUTDIR, sample=SAMPLE_SET, num=SET_NUMS),
    # raw_zip = expand("{outputdir}00_fastqc/{sample}_L001_{num}_001_fastqc.zip", outputdir=OUTPUTDIR, sample=SAMPLE_SET, num=SET_NUMS),
    # raw_multi_html = OUTPUTDIR + "00_fastqc/raw_multiqc.html",    
    # # reference for cellranger
    # out_ref = REF + config["reference"]["ref_cellranger"] + "/fasta/genome.fa",
    # # Cellranger count
    # out_cellranger = OUTPUTDIR + "00_cellranger/" + config["sample"]["sname"] + "/outs/web_summary.html",
    seurat_report = OUTPUTDIR + "01_seurat/seurat.html",

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: ENVDIR + "fastqc.smk"
include: ENVDIR + "cellranger.smk"
include: ENVDIR + "seurat.smk"
