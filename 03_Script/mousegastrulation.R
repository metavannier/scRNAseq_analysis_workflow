#..............................................................
# Preparation of the folowing reference for sims :
# Mouse gastrulation reference :
# doi:10.1038/s41586-019-0933-9
#..............................................................

#-------------------------------------
# Library
#-------------------------------------

library(magick)
library(BiocManager)
library(reshape2)
library(Matrix)
library(data.table)
library(scran)
library(igraph)
library(scuttle)
library(scater)

BiocManager::install("MouseGastrulationData")
library(MouseGastrulationData)

#-------------------------------------
# Path to files / folders
#-------------------------------------

DIRECTORY = getwd()
OUTPUTDIR = file.path((DIRECTORY), "05_Output")

SAMPLE_ID = snakemake@params[["sample_id"]]

MOUSEGASTRULATION_SAMPLES =  snakemake@params[["mousegastrulation_samples"]]

OUTPUT_NAME_REF_METADATA = snakemake@params[["output_name_ref_metadata"]]
OUTPUT_NAME_REF_MATRIX = snakemake@params[["output_name_ref_matrix"]]
OUTPUT_NAME_MATRIX = snakemake@params[["output_name_matrix"]]
TEXT_OUTPUT = snakemake@output[["data_for_sims_output"]]

STEP2 = "02_seurat/"
STEP3 = "03_sims"

#-------------------------------------
# Path to files / folders
#-------------------------------------

##  Loading the samples from the atlas 
MOUSEGASTRULATION_SAMPLES <- type.convert(MOUSEGASTRULATION_SAMPLES, dec=".")

sce <- EmbryoAtlasData(samples = c(MOUSEGASTRULATION_SAMPLES))

#-------------------------------------
# Normalisation : we follow the process from Pijuan-Sala et al. 2019 
# (doi:10.1038/s41586-019-0933-9)
#-------------------------------------


#### DANS LE GIT HUB, VOIR SI A FAIRE SACHANT QUE LES SIZEFACTOR SONT DANS L'OBJET sce
head(sizeFactors(sce))
## Size factors for normalisation with scran.
sce = sce[calcAverage(sce)>0.1,]

# For pre-clustering, we use scran's `quickCluster` function, using the `igraph` method. We specify a maximum cluster size of 3000 cells and a minimum cluster size of 100 cells.
clusts = as.numeric(quickCluster(sce, method = "igraph", min.size = 100))
#now run the normalisation
#number of cells in each cluster should be at least twice that of the largest 'sizes'
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
compute_sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)
### Pour voir si on trouve la même normalisation que dans l'objet :
head(compute_sce)
#####

# Normalisation using use scuttle to get lognormcounts function
sce <- logNormCounts(sce)
normcounts_atlas <- normcounts(sce)

#-------------------------------------
# Metadata reference file
#-------------------------------------

## Produce the reference metadata matrice
reference_metadata <- cbind(sce$cell,sce$celltype)

#-------------------------------------
# Wtrite the two new matrix +
# create directory
#-------------------------------------

dir.create(file.path(OUTPUTDIR, STEP3, SAMPLE_ID))

class(sce) <- "numeric"
write.csv(normcounts_atlas, file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_MATRIX,".csv")))

fwrite(x = reference_metadata, file = file.path(OUTPUTDIR, STEP3, SAMPLE_ID, paste0(OUTPUT_NAME_REF_METADATA,".csv")))

#-------------------------------------
# create the output file for the snakemake rule
#-------------------------------------

output_file<-file(TEXT_OUTPUT)
writeLines(c("Files preparation for the Mouse Gastrulation reference for SIMS finished (CSV format)"), output_file)
close(output_file)
