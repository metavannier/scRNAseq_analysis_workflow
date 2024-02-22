This html page present standard pre-processing workflow by DropletUtils, Seurat and scran.

DropletUtils
- Swapped molecules were excluded using the swappedDrops() function from DropletUtils (Griffiths et al. 2018).
- Cell-containing droplets were called using the emptyDrops() function from DropletUtils (Lun et al. 2019).

Seurat:
- QC and selecting cells for further analysis
- Normalizing the data
- Identification of highly variable features
- Scaling the data
- Perform linear dimensional reduction
- Determine the ‘dimensionality’ of the dataset
- Cluster the cells
- Assigning cell type identity to clusters

Scran:
- Size factors were computed using the computeSumFactors() function from scran (Lun, Bach, and Marioni 2016).
- Putative doublets were identified and excluded using the doubletCells() function from scran.
- Cytoplasm-stripped nuclei were also excluded.
- Batch correction was performed in the principal component space with fastMNN() from scran (Haghverdi et al. 2018).