# ----------------------------------------------
# This script is used to launch the compilation
# of the RMarkdown report independently 
# of Rstudio interface
# ----------------------------------------------

WORKING_DIR = getwd()

STEP1 =  "02_seurat"
STEP2 = "rm_doublet"

SCRIPTDIR = file.path( WORKING_DIR, "03_Script")
OUTPUTDIR = file.path( WORKING_DIR, "05_Output")

# General
RULES_OUTPUT <- snakemake@output[["doublets_output"]]
SEURAT_DATA_PATH_LIST <- snakemake@params[["seurat_data"]]
REPORT_LIST <- snakemake@output[["doublets_report"]]
SAMPLE_ID_LIST <- snakemake@params[["sample_id"]]

len <- length(SAMPLE_ID_LIST)
i = 1

while(i <= len)
{
    
    SEURAT_DATA_PATH <- SEURAT_DATA_PATH_LIST[i]
    REPORT <- REPORT_LIST[i]
    SAMPLE_ID <- SAMPLE_ID_LIST[i]

    rmarkdown::render( input = file.path(SCRIPTDIR, "rm_doublets.Rmd"),
                    output_dir = file.path(OUTPUTDIR, STEP2, SAMPLE_ID),
                    output_file = REPORT,
                    quiet = FALSE)
    i = i + 1
}

# output_file<-file(RULES_OUTPUT)
# writeLines(c("Rules rm_doublet FINISHED"), output_file)
# close(output_file)

    