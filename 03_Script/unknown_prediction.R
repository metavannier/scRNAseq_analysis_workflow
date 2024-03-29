#-------------------------------------
# This script say if a cell is considered as unknown
# depending of the prediction score of sims
#-------------------------------------

library(data.table)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")

SAMPLE_ID = snakemake@params[["sample_id"]]
MATRIX = snakemake@params[["matrix"]]
THRESHOLD = snakemake@params[["threshold"]]
PRED_FILTERED = snakemake@params[["pred_filtered"]]

STEP3 = "03_sims/"

# Load matrix with sims prediction 
pred <- fread(file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, MATRIX,"_prediction.csv")))

# Assigning Cells, Predictions and probabilities to finalpred
finalpred <- pred[, c("Cells", "first_pred", "second_pred","first_prob", "second_prob"), drop = FALSE]

# Assigning predictions and probabilities to finalpred
diff_prob <- pred$first_prob - pred$second_prob

# Replace the value of first_pred if diff_prob is less than the threshold
finalpred$first_pred[diff_prob < THRESHOLD] <- "unknown"

#-------------------------------------
# create the predicted matrix with the
# output file for the snakemake rule
#-------------------------------------

write.table(file=file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_",PRED_FILTERED)), finalpred ,sep=",", quote=F, row.names=F, col.names=T)

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

TEXT_OUTPUT <- snakemake@output[["unknown_prediction_output"]]

output_file<-file(TEXT_OUTPUT)
writeLines(c("Rules unknown prediction finished"), output_file)
close(output_file)
