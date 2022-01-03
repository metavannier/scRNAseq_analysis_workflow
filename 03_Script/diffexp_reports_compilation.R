# ##########################################################################
# This script is used to launch the compilation of the RMarkdown report
# independently of Rstudio interface
# ##########################################################################

STEP = "03_diffexp"
WORKING_DIR = getwd()
SCRIPTDIR = file.path( WORKING_DIR, "03_Script")
OUTPUTDIR = file.path( WORKING_DIR, "05_Output")
REPORT = snakemake@output[["diffexp_report"]]
NPROJ = snakemake@params[["nproj"]]

OUTPUTFILE = paste(NPROJ, "diffexp_report.html", sep = "_")

rmarkdown::render( input = file.path(SCRIPTDIR, "diffexp.Rmd"),
                   output_dir = file.path(OUTPUTDIR, STEP),
                   output_file  = OUTPUTFILE,
                   quiet = FALSE)
