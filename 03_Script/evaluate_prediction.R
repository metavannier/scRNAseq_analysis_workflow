#-------------------------------------
# This script evaluate the prediction obtain
# by sims after the unknown assignation of
# cells with low probability score
#-------------------------------------

library(data.table)
# library(SingleCellExperiment)
# BiocManager::install("MouseGastrulationData")
# library(MouseGastrulationData)
library(ggplot2)

#-------------------------------------
# Path to files / folders
#-------------------------------------
DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")
REF = file.path((DIRECTORY), "01_Reference")

SAMPLE_ID = snakemake@params[["sample_id"]]
MATRIX = snakemake@params[["matrix"]]
PRED_FILTERED = snakemake@params[["pred_filtered"]]
MOUSEGASTRULATION_SAMPLES =  snakemake@params[["mousegastrulation_samples"]]
CONFUSION_MATRIX =  snakemake@params[["confusion_matrix"]]

STEP3 = "03_sims/"

#-------------------------------------
# Path to files / folders
#-------------------------------------

# Load matrix with sims prediction and unknown cells
cell_predict <- fread(file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_",PRED_FILTERED)))

# Load the atlas with prediction and unknown cells
# MOUSEGASTRULATION_SAMPLES <- type.convert(MOUSEGASTRULATION_SAMPLES, dec=".", as.is = TRUE)
# sce <- EmbryoAtlasData(samples = c(MOUSEGASTRULATION_SAMPLES))
# cell_atlas <- DataFrame(Cells=colData(sce)$cell,CellType_atlas=colData(sce)$celltype)

# write.table(file=file.path("01_Reference/atlas/mousegastrulation_atlas_cell_annotation"), cell_atlas ,sep=",", quote=F, row.names=F, col.names=T)

# Create a new dataframe by merging based on the 'Cells' column
cell_atlas <- fread(file = file.path(REF, "atlas", paste0(SAMPLE_ID, "_atlas_cell_annotation")))

# Replace NA values with "unknown" in the "first_pred" column
cell_atlas$CellType_atlas <- replace(cell_atlas$CellType_atlas, is.na(cell_atlas$CellType_atlas), "unknown")

cell_eval <- merge(cell_predict, cell_atlas, by = "Cells")

# Assigning Cells, Predictions and probabilities to finalpred
cell_eval_matrix <- cell_eval[, c("Cells", "first_pred", "CellType_atlas"), drop = FALSE]

#-------------------------------------
# Create a confusion matrix
#-------------------------------------

confusion_matrix <- table(cell_eval_matrix$first_pred, cell_eval_matrix$CellType_atlas)

# Convert the confusion matrix table to a data frame
confusion_df <- as.data.frame.table(confusion_matrix)

# Subset the confusion_df to rows where Var2 is "unknown" and sum the counts
unknown_sum <- sum(confusion_df[confusion_df$Var2 == "unknown", "Freq"])
print(unknown_sum)

# Compute total number of cells
total_sum <- sum(confusion_df$Freq)
print(total_sum)

# Remove "unknown" row and column from the confusion matrix
confusion_matrix_filtered <- confusion_matrix[-which(rownames(confusion_matrix) == "unknown"), -which(colnames(confusion_matrix) == "unknown")]

# Calculate scores
accuracy <- sum(diag(confusion_matrix_filtered)) / sum(confusion_matrix_filtered)
precision <- diag(confusion_matrix_filtered) / rowSums(confusion_matrix_filtered)
recall <- diag(confusion_matrix_filtered) / colSums(confusion_matrix_filtered)
f1_score <- 2 * (precision * recall) / (precision + recall)

# Calculate the mean of precision, recall, and F1 score
mean_precision <- mean(precision, na.rm = TRUE)
mean_recall <- mean(recall, na.rm = TRUE)
mean_f1_score <- mean(f1_score, na.rm = TRUE)

# Create a data frame for scores
score_df <- data.frame(metric = c("Accuracy", "Mean Precision", "Mean Recall", "Mean F1 Score"),
                       value = c(accuracy, mean_precision, mean_recall, mean_f1_score))

print(score_df)

# Create the plot
confusion_plot <- ggplot(data = confusion_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "Predicted", y = "Atlas", title = "Confusion Matrix") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

# Print the plot
ggsave(file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_",CONFUSION_MATRIX)), confusion_plot, width = 20, height = 18)

# -------------------------------------
# create the output file for the snakemake rule
# -------------------------------------

# TEXT_OUTPUT <- snakemake@output[["evaluate_prediction_output"]]

# output_file<-file(TEXT_OUTPUT)
# writeLines(c("Rules evaluate prediction finished"), output_file)
# close(output_file)

