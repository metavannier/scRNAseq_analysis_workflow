#..............................................................
# merge_rds : Merging the count of samples from rds files 
#..............................................................

#-------------------------------------
# Path to files / folders
#-------------------------------------

DIRECTORY = getwd()
OUTPUTDIR = "05_Output"
SAMPLE_ID = snakemake@params[["sample_id"]]
STEP2 = "02_seurat/"
STEP = "04_velocity/"
TEXT_OUTPUT = snakemake@output[["merge_rds_output"]]

#-------------------------------------
# Loading library
#-------------------------------------
library(Seurat)

## Data preparation 
# Merge the rds of each samples after seurat filtration

# Read the RDS files into a list
sep_data.list <- lapply(file.path(OUTPUTDIR,STEP2, SAMPLE_ID, paste0(SAMPLE_ID,"_filtered_seurat_object.rds")), readRDS)

# Merge rds
# SO <- merge(x = sep_data.list[[1]], y = sep_data.list[2:length(sep_data.list)], merge.data = T)

# Merge the data frames in the list
merged_data <- Reduce(function(x, y) merge(x, y), sep_data.list)

# metadata <- merged_data@meta.data
# unique_orig_idents <- unique(metadata$orig.ident)
# print(unique_orig_idents)
## Remove/show metadata (To do?)
merged_data$seurat_clusters <- NULL
merged_data@meta.data[grep("pANN", colnames(merged_data@meta.data), value = TRUE)] <- NULL
merged_data@meta.data[grep("RNA_snn", colnames(merged_data@meta.data), value=TRUE)] <- NULL
# head(merged_data@meta.data)
# Save merge
saveRDS(merged_data, file = file.path( OUTPUTDIR, STEP, "merged_seurat_object.rds"))

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

output_file<-file(TEXT_OUTPUT)
writeLines(c("Merging step finished"), output_file)
close(output_file)

sessionInfo(package = NULL)