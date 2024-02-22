# ########################################################
# This script exclude the swapped molecules  and 
# cell-containing droplets with DropletUtils
# ########################################################

## @knitr swappedDrops

#........................................
# Debug 
#........................................

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;

#........................................
# Path to files / folders
#........................................


DIRECTORY = getwd()
OUTPUTDIR = file.path(dirname(DIRECTORY), "05_Output")
STEP1 = "01_cellranger/"
STEP2 = "02_seurat/"

source(file.path(dirname(DIRECTORY), "03_Script/00_general_deps.R"))

in_data_dir = file.path( OUTPUTDIR, STEP1)
samples <- dir(in_data_dir)

ncores = 3
mcparam = MulticoreParam(workers = ncores)
register(mcparam)

#........................................
# Read data
#........................................

mol_loc = paste0(in_data_dir,SAMPLE_ID_LIST,"/count/sample_molecule_info.h5")
mtx_loc = paste0(in_data_dir,SAMPLE_ID,CELL_RANGER_COUNT_PATH, "matrix_unswapped.mtx")
bc_loc = paste0(in_data_dir,SAMPLE_ID,CELL_RANGER_COUNT_PATH, "barcodes_unswapped.tsv")
gene_loc = paste0(in_data_dir,SAMPLE_ID,CELL_RANGER_COUNT_PATH, "genes_unswapped.tsv")

unswapped = swappedDrops(mol_loc, get.swapped = TRUE)

ratios = sapply(1:length(unswapped$cleaned), function(i){
  sum(unswapped$swapped[[i]])/(sum(unswapped$cleaned[[i]]) + sum(unswapped$swapped[[i]]))
})

names(ratios) = SAMPLE_ID_LIST

mat = writeMM(unswapped$cleaned[[as.numeric(INCREMENT_SAMPLE)]], file = mtx_loc)

write.table(colnames(unswapped$cleaned[[as.numeric(INCREMENT_SAMPLE)]]), file = bc_loc, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(rownames(unswapped$cleaned[[as.numeric(INCREMENT_SAMPLE)]]), file = gene_loc, col.names = FALSE, row.names = FALSE, quote = FALSE)

ggplot(data.frame(ratios = ratios, sample = SAMPLE_ID), aes(x = sample, y = ratios)) +
  geom_bar(stat = "identity", fill = "grey20") +
  labs(x = "Sample", y = "Fraction of molecules excluded")

## @knitr emptyDrops

#........................................
# Cell calling
#........................................

# Load the data
###

DIRECTORY = getwd()
OUTPUTDIR = file.path(dirname(DIRECTORY), "05_Output")
STEP1 = "01_cellranger/"
STEP2 = "02_seurat/"
source(file.path(dirname(DIRECTORY), "03_Script/00_general_deps.R"))
in_data_dir = file.path( OUTPUTDIR, STEP1)
samples <- dir(in_data_dir)
mol_loc = paste0(in_data_dir,SAMPLE_ID,"/count/sample_molecule_info.h5")
# mtx_loc = paste0(in_data_dir,SAMPLE_ID,CELL_RANGER_COUNT_PATH, "matrix_unswapped.mtx")
bc_loc = paste0(in_data_dir,SAMPLE_ID,CELL_RANGER_COUNT_PATH, "barcodes_unswapped.tsv")
gene_loc = paste0(in_data_dir,SAMPLE_ID,CELL_RANGER_COUNT_PATH, "genes_unswapped.tsv")
####

mtx_loc = paste0(in_data_dir,SAMPLE_ID_LIST,CELL_RANGER_COUNT_PATH, "matrix_unswapped.mtx")

matrices = bplapply(mtx_loc, readMM)

## Seems doe'nt work ##
mat = matrices[[INCREMENT_SAMPLE]]
bcs = read.table(bc_loc, header = FALSE, stringsAsFactors = FALSE)[,1]
# #correct barcode sample number
bcs = paste0(bcs, "-", i)

targets = mat[, Matrix::colSums(mat)!=0 & Matrix::colSums(mat) < 100]
sums = Matrix::colSums(targets)
reads = as.matrix(table(sums) * as.numeric(names(table(sums))))

ggplot(data.frame(lib = reads, n = as.numeric(rownames(reads))), aes(x = n, y= lib)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Barcode library UMI count", y = "Total contribution to background vector") +
  scale_x_continuous(breaks = seq(0,100, 5))

bg_vec = Matrix::rowSums(targets)

####

set.seed(42)

mtx_loc = paste0(in_data_dir,SAMPLE_ID_LIST,CELL_RANGER_COUNT_PATH, "matrix_unswapped.mtx")
matrices = bplapply(mtx_loc, readMM)

outs = lapply(matrices, emptyDrops, niters = 20000, ignore = 4999, BPPARAM = mcparam, lower = 100, retain = Inf)

sessionInfo()
