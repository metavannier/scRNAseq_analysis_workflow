----------------------------------
scRNA-seq analysis of gastruloids.
----------------------------------

This workflow performs a Single-Cell RNA-seq analysis of gastruloids at 4 different time points.

---------------------
Materials and methods
---------------------

- The quality of the raw reads are assessed using `FastQC v0.12.1 toolkit <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ (Andrews et al., 2010). Nextera Transposase Sequence are removed and low quality reads are trimmed using Trimmomatic v0.39.

- Raw files were processed with CellRanger v7.1.0 (12/07/2022) using default mapping arguments, with reads mapped to the mouse GRCm38 primary sequence assembly and counted with the GRCm38.98 optimized mouse annotation v2.0 from the `ReferenceEnhancer R package <https://doi.org/10.1038/s41592-023-02003-w>`_ (Pool et al. Nature methods 2023). Cell Ranger multi pipeline was run with the intron mode used to generate exonic and intronic read mapping based gene-cell matrices.

- QC steps were performed with `Seurat version 4.0.3 <https://doi.org/10.1038/nbt.4096>`_ (Butler et al. Nature Biotechnology 2018). loaded into R version 4.3.2. Low quality cells were removed with the following thresholds: Unique feature counts (UFC) > 1700 and < 7900 & mitochondrial counts (MC) < 10% & Total RNA counts (TRC) < 55000 for the 72h dataset, UFC > 1700 and < 7500 & MC < 10% & TRC < 42000 for the 80h dataset, UFC > 1400 and < 7700 & MC < 8% & TRC < 42000 for the 86h dataset, UFC > 1400 and < 7200 & MC < 8% & TRC < 40000 for the 96h dataset.

- We normalized the feature expression measurments for each cell with the same method as the paper `mouse atlas <https://doi.org/10.1038/s41586-019-0933-9>`_ (Andrews, 2010) :

Transcriptome size factors were calculated for each dataset separately, using ‘computeSumFactors’ from the scran R package. 
Cells were pre-clustered with the ‘quickCluster’ function using the parameter ‘method=igraph’ (using the scran R package), and minimum and maximum cluster sizes of 100 and 3,000 cells, respectively. 
Raw counts for each cell were divided by their size factors, and the resulting normalized counts were used for further processing.

- We use the `mouse gastrulation and early organogenesis atlas <https://doi.org/10.1038/s41586-019-0933-9>`_ (Pijuan-Sala et al.) without the non annotated cells as reference with the software `SIMS <https://doi.org/10.1101%2F2023.02.28.529615>`_ (Gonzalez-Ferrer et al. 2023). The training step was stoped after 34 epoch by visualysing the loss and acccuracy curve evolution. We define as "Unknown" a cell with a difference less than 0.10 between the probability measure of the first and the second classification.

Optimizing the transcriptomic Reference
=======================================

We use the `ReferenceEnhancer R package <https://doi.org/10.1038/s41592-023-02003-w>`_ to recovered missing gene expression data by optimizing the reference transcriptome for scRNA-seq through recovering :

A) Poor annotation of 3' gene ends.

B) Issues with intronic read incorporation.
  
C) Gene overlap-derived read loss.

- "We demonstrate, with a diverse collection of mouse and human tissue data, that reference optimization can substantially improve cellular profiling resolution and reveal missing cell types and marker genes."

Normalization
=============

The size-factor normalization reduces the technical variation from sequencing depth. 
The size factor for each cell represents the estimate of the relative bias in that cell, so division of its counts by its size factor should remove that bias.
Normalization methods used in bulk RNA-seq data (e.g. DESeq2 size-factors) are not appropriate for sc-RNAseq data due to the large numbers of zeros in the data.

We use the same method as the `Pijuan-Sala et al. paper <https://doi.org/10.1038/s41586-019-0933-9>`_  consisting on pooling similar cells together, one can obtain enough data to estimate a size factor :

- Clusters cells to group cells with similar gene expression.

- Pool cells within these clusters and calculate size factors.

- Repeat step 2 with varying sets of cells and pool sizes.

- Derive cell-specific size factors using linear algebra.

Cell classification
===================

We use the software `SIMS <https://doi.org/10.1101%2F2023.02.28.529615>`_ for the Cell classification.

SIMS :  A Deep Learning Label Transfer Tool for Single-Cell RNA Sequencing Analysis

- Tools based on `TabNet <https://arxiv.org/abs/1908.07442>`_ as the classifier component. Optimized for tabular data, can work well with over 45K features. TabNet is a transformer-based neural network but don't need pre-training.

- SIMS model takes advantage of lack of expression, and fluctuations of expression levels of the whole transcriptome to learn and identify cell labels.

- Temperature Scaling, a post-processing step of the train network that provides the users with a calibrated probability measure for the classification of each cell.