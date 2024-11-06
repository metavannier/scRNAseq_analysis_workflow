# ----------------------------------------------
# This script is used to launch the compilation
# of the RMarkdown report for the velocity analyses with URD 
# independently of Rstudio interface
# ----------------------------------------------

WORKING_DIR = getwd()

STEP = "04_velocity"

SCRIPTDIR = file.path( WORKING_DIR, "03_Script")
OUTPUTDIR = file.path( WORKING_DIR, "05_Output")


rmarkdown::render( input = file.path(SCRIPTDIR, "velocity/urd.Rmd"),
                output_dir = file.path(OUTPUTDIR, STEP),
                output_file = "urd_report.html",
                quiet = FALSE)


output_file<-file(snakemake@output[["urd_output"]])
writeLines(c("Rule urd FINISHED"), output_file)
close(output_file)