#!/usr/bin/env Rscript 

library(devtools)
install_github("PoolLab/ReferenceEnhancer", ref = "1b9d1cb")

library(ReferenceEnhancer)

args <- commandArgs(trailingOnly = TRUE)

ref_gtf <- args[1]

genome_annotation <- LoadGtf(ref_gtf)

# Identify all overlapping genes based on the ENSEMBL/10x Genomics genome annotation file (GTF), rank-order them according to the number of gene overlaps.
gene_overlaps <- IdentifyOverlappers(genome_annotation)

print("OverlapResolutions")
# Generate recommended actions for overlapping genes based on original genome annotation .gtf file and a list of overlapping genes.
# Does not work : problem with first column?
OverlapResolutions(genome_annotation, gene_overlaps, c("Rik$", "^Gm"))

print("IsolateIntergenicReads")
# Extract intergenic reads from Cell Ranger aligned bam file.
IsolateIntergenicReads("S004647_to_S004650/outs/er_sample_outs/S004648/count/sample_alignments.bam", "S004647_to_S004650/outs/er_sample_outs/S004648/count/sample_alignments.bam.bai")

print("GenerateGeneLocationBed")
GenerateGeneLocationBed(genome_annotation)

GenerateExtensionCandidates() 

OptimizedAnnotationAssembler(ref_gtf, "overlapping_gene_list.csv", "gene_extension_candidates.csv", "rename_genes.csv")

sessionInfo()