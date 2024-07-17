#-------------------------------------
# This script say if a cell is considered as unknown
# depending of the prediction score of sims
#-------------------------------------

library(data.table)
library(dplyr)

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

# Compute the differences
differences <- pred$first_prob - pred$second_prob

# Perform the Shapiro-Wilk test for normality
shapiro_test_result <- shapiro.test(differences)
shapiro_p_value <- shapiro_test_result$p.value

# Perform the Kolmogorov-Smirnov test for normality
ks_test_result <- ks.test(differences, "pnorm", mean=mean(differences), sd=sd(differences))
ks_p_value <- ks_test_result$p.value

# Interpretation and testing
alpha <- 0.05

if (shapiro_p_value > alpha & ks_p_value > alpha) {
  # Perform a one-sample t-test
  t_test_result <- t.test(differences, mu=0, alternative="greater")
  # Print the results of the t-test
  print(t_test_result)
} else { 
  # Perform the Wilcoxon signed-rank test
  wilcox_test_result <- wilcox.test(differences, mu=0, alternative="greater")
  # Print the results of the Wilcoxon signed-rank test
  print(wilcox_test_result)
}

png(filename = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, MATRIX, "_score_differences_histogram.png")))
hist(differences, breaks=10, main="Distribution of Score Differences", xlab="Difference (Score1 - Score2)", col="skyblue", border="black")
dev.off()

# Create a Q-Q plot
qqnorm(differences)
qqline(differences, col="red")

# Save the Q-Q plot to a file
png(filename = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, MATRIX, "qq_plot.png")))
qqnorm(differences)
qqline(differences, col="red")
dev.off()

# Compute the confidence interval for the differences
# Assuming a 95% confidence interval
alpha <- 0.05
n <- length(differences)
mean_diff <- mean(differences)
std_error <- sd(differences) / sqrt(n)
z_score <- qnorm(1 - alpha / 2)  # Z-score for 95% confidence interval

# Calculate the lower bound of the confidence interval
lower_bound <- mean_diff - z_score * std_error
print(lower_bound)
# Apply the threshold to update predictions
finalpred$final_pred <- ifelse(differences <= lower_bound, "unknown", finalpred$first_pred)

#-------------------------------------
# create the predicted matrix with the
# output file for the snakemake rule
#-------------------------------------

write.table(file=file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_",PRED_FILTERED)), finalpred ,sep=",", quote=F, row.names=F, col.names=T)

#-------------------------------------
# Create a table to record the counts 
# of the two first an second prediction 
# when unknown
#-------------------------------------

# Create a new column with sorted combinations of first_pred and second_pred
finalpred_unknown <- finalpred %>%
  filter(final_pred == "unknown") %>%
  mutate(binome = ifelse(first_pred < second_pred,
                         paste(first_pred, second_pred, sep = "-"),
                         paste(second_pred, first_pred, sep = "-")))

# Count occurrences of each binome
binome_table <- finalpred_unknown %>%
  count(binome) %>%
  rename(binome_count = n) %>%
  filter(binome_count > 10) %>%
  arrange(desc(binome_count))

write.table(file=file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(SAMPLE_ID, "_binome_unknown.csv")), binome_table ,sep=",", quote=F, row.names=F, col.names=T)


#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

TEXT_OUTPUT <- snakemake@output[["unknown_prediction_output"]]

output_file<-file(TEXT_OUTPUT)
writeLines(c("Rules unknown prediction finished"), output_file)
close(output_file)
