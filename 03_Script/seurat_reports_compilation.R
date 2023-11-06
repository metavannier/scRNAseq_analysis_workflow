# ----------------------------------------------
# This script is used to launch thr compilation
# of the RMarkdown report independently 
# of Rstudio interface
# ----------------------------------------------

WORKING_DIR = getwd()

STEP = "02_seurat"

SCRIPTDIR = file.path( WORKING_DIR, "03_Script")
OUTPUTDIR = file.path( WORKING_DIR, "05_Output")

# General
SEURAT_OUTPUT <- snakemake@output[["seurat_output"]]
SC_DATA_PATH_LIST <- snakemake@params[["sc_data"]]
REPORT_LIST <- snakemake@output[["seurat_report"]]
SAMPLE_ID_LIST <- snakemake@params[["sample_id"]]
PLOT_RASTER_NBCELLS_THRESHOLD <- snakemake@params[["plot_raster_nbcells_threshold"]]
CELL_RANGER_COUNT_PATH <- snakemake@params[["cell_ranger_count_path"]]
RUN_DEMULTIPLEX <- snakemake@params[["run_demultiplex"]]

# Filtering
QC_EXPLORATION_MODE <- snakemake@params[["qc_exploration_mode"]]
MIN_CELLS <-snakemake@params[["min_cells"]]
MIN_FEATURES <- snakemake@params[["min_features"]]
FILTER_UMI_MIN_LIST <- snakemake@params[["filter_umi_min"]]
FILTER_UMI_MAX_LIST <- snakemake@params[["filter_umi_max"]]
FILTER_FEATURE_MIN_LIST <- snakemake@params[["filter_feature_min"]]
FILTER_FEATURE_MAX_LIST <- snakemake@params[["filter_feature_max"]]
FILTER_PERCENT_MT_LIST <- snakemake@params[["filter_percent_mt"]]
FILTER_PERCENT_MT_MIN_LIST <- snakemake@params[["filter_percent_mt_min"]]
FILTER_PERCENT_RB_LIST <- snakemake@params[["filter_percent_rb"]]
PATTERN_MT <- snakemake@params[["pattern_mt"]]
PATTERN_RB <- snakemake@params[["pattern_rb"]]

# Normalization
NORM_METHOD <- snakemake@params[["norm_method"]]
NORM_SCALE_FACTOR <- snakemake@params[["norm_scale_factor"]]

FEATURE_SELECT_METHOD <- snakemake@params[["feature_select_method"]]
VARIABLE_FEATURES <- snakemake@params[["variable_features"]]
VARIABLE_FEATURES_SHOWTOP <- snakemake@params[["variable_features_showtop"]]

# PCA parameters
DIMS <- snakemake@params[["dims"]]
PCA_NPC <- snakemake@params[["pca_npc"]]
PCA_PLOTS_NBDIMS <- snakemake@params[["pca_plots_nbdims"]]
PCA_PLOTS_NBFEATURES <- snakemake@params[["pca_plot_nbfeatures"]]

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS <- snakemake@params[["dimreduc_use_pca_nbdims"]]

# Cluster identification parameters
FINDCLUSTERS_USE_PCA_NBDIMS <- snakemake@params[["findclusters_use_pca_nbdims"]]
FINDNEIGHBORS_K <- snakemake@params[["findneighbors_k"]]
FINDCLUSTERS_RESOLUTION <- snakemake@params[["findclusters_resolution"]]
FINDCLUSTERS_ALGORITHM <- snakemake@params[["findclusters_algorithm"]]


len <- length(SAMPLE_ID_LIST)
i = 1

while(i <= len)
{
    
    SC_DATA_PATH <- SC_DATA_PATH_LIST[i]
    REPORT <- REPORT_LIST[i]
    SAMPLE_ID <- SAMPLE_ID_LIST[i]
    FILTER_UMI_MIN <- FILTER_UMI_MIN_LIST[i]
    FILTER_UMI_MAX <- FILTER_UMI_MAX_LIST[i]
    FILTER_FEATURE_MIN <- FILTER_FEATURE_MIN_LIST[i]
    FILTER_FEATURE_MAX <- FILTER_FEATURE_MAX_LIST[i]
    FILTER_PERCENT_MT <- FILTER_PERCENT_MT_LIST[i]
    FILTER_PERCENT_MT_MIN <- FILTER_PERCENT_MT_MIN_LIST[i]
    FILTER_PERCENT_RB <- FILTER_PERCENT_RB_LIST[i]

    rmarkdown::render( input = file.path(SCRIPTDIR, "seurat.Rmd"),
                    output_dir = file.path(OUTPUTDIR, STEP, SAMPLE_ID),
                    output_file = REPORT,
                    quiet = FALSE)
    i = i + 1
}

output_file<-file(SEURAT_OUTPUT)
writeLines(c("Rule seurat FINISHED"), output_file)
close(output_file)

    