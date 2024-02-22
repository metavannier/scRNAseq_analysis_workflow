## Removes doublets

## @knitr computescores

TEXT_OUTPUT = snakemake@output[["doublets_output"]]
SEURAT_DATA = snakemake@params[["seurat_data"]]
NORMALISATION = snakemake@params[["normalisation"]]
# Computing doublet scores

# sce <- readRDS(file.path( OUTPUTDIR, STEP1, SAMPLE_ID, paste0(SAMPLE_ID, "_filtered", NORMALISATION, "_seurat_object.rds")))
# # sce = fread(file.path( OUTPUTDIR, STEP1, SAMPLE_ID, paste0(SAMPLE_ID, "_", NORMALISATION, "_normalized_matrix.csv")))
counts <- read.table(file.path( OUTPUTDIR, STEP1, SAMPLE_ID, paste0(SAMPLE_ID, "_", NORMALISATION, "_normalized_matrix.csv")), header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "")
genes = counts[,1]
cell = colnames(counts)


counts = as.matrix(as.numeric(counts))

rownames(counts) = genes
colnames(counts) = cell
counts=counts[-c(1),-c(1)]
head(counts)

sce = SingleCellExperiment(assays = list("counts" = counts))
# Convert to SingleCellExperiment
# sce <- as.SingleCellExperiment(sce)

require(biomaRt)
trend = scran::trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
# decomp = scran::decomposeVar(sce, fit = trend)
# decomp = decomp[decomp$mean > 1e-3,]
# xist = "ENSMUSG00000086503"
# mouse_ensembl = useMart("ensembl")
# mouse_ensembl = useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)
# gene_map = getBM(attributes=c("ensembl_gene_id", "chromosome_name"), filters = "ensembl_gene_id", values = rownames(decomp), mart = mouse_ensembl)
# ychr = gene_map[gene_map[,2] == "Y", 1]
# other = c("tomato-td") #for the chimera
# decomp = decomp[!rownames(decomp) %in% c(xist, ychr, other),]
# decomp$FDR = p.adjust(decomp$p.value, method = "fdr")






# hvg.list = getHVGs(sce)

# names(hvg.list) = unique(meta$sample)
# set.seed(42)

# scores_hvgs = doubletCells(sce, approximate = TRUE, subset.row = rownames(sce) %in% hvg.list)

# scores_hvgs = do.call(c, scores_hvgs)

# scores = scores_hvgs

# print(scores)

## @knitr info
sessionInfo()

# create the output file for the snakemake rule

output_file<-file(TEXT_OUTPUT)
writeLines(c("Removes doublets step finished"), output_file)
close(output_file)