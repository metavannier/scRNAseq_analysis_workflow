# ##########################################################################
# This script is used to launch the compilation of the RMarkdown report
# independently of Rstudio interface
# ##########################################################################

WORKING_DIR = getwd()

STEP = "02_seurat"

SCRIPTDIR = file.path( WORKING_DIR, "03_Script")
OUTPUTDIR = file.path( WORKING_DIR, "05_Output")
NPROJ = snakemake@params[["nproj"]]
REPORT = snakemake@output[["seurat_report"]]
SC_DATA_PATH = snakemake@input[["sc_data"]]
OUTPUTFILE = paste(NPROJ, "seurat_report.html", sep = "_")

rmarkdown::render( input = file.path(SCRIPTDIR, "seurat.Rmd"),
                   output_dir = file.path(OUTPUTDIR, STEP),
                   output_file  = OUTPUTFILE,
                   quiet = FALSE)
