# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

# ..........................................................................................................
## @knitr heterogeneity_dimReduc
# ..........................................................................................................

#........................................
# Dimensional reduction (tSNE / UMAP)
#........................................
# Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
seurat_obj = RunUMAP( seurat_obj, dims = 1:pcs);

# Run t-SNE dimensionality reduction on selected features
seurat_obj = RunTSNE( seurat_obj, dims = 1:pcs);

# Save resulting coordinates for all cells as 'tsv' files
# write.table( Embeddings(seurat_obj, reduction = "umap"), 
#              file= file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_cellsCoordinates_umap.tsv")), 
#              quote = FALSE, 
#              row.names = TRUE, 
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");

# write.table( Embeddings(seurat_obj, reduction = "tsne"), 
#              file= file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "cellsCoordinates_tsne.tsv")), 
#              quote = FALSE, 
#              row.names = TRUE, 
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");

# ..........................................................................................................
## @knitr dimreduc_ggplot_covariables
# ..........................................................................................................

# Plot the cells percentage of ribosomal genes in UMAP
FeaturePlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
             features = "percent.rb") +
              scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
              ggtitle( "Map of cells with level of\n percentage of ribosomal genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by %of ribosomal genes or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                                  group.by = "outlier.percent.rb") + 
                 ggtitle( "Map of cells by filter status\n on %ribosomal gene")
      )
}

# Plot the cells percentage of mitochondrial genes in UMAP  
FeaturePlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "percent.mt") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with level of percentage of\n mitochondrial genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by %of mitochondrial genes or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.percent.mt") + 
           ggtitle( "Map of cells by filter status\n on %mitochondrial gene")
  )
}

# Plot the cells RNA counts in UMAP
FeaturePlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nCount_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with RNA counts")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by number of UMI count or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.nCount_RNA") + 
           ggtitle( "Map of cells by filter status \non nCount_RNA value")
  )
}

# Plot the cells genes counts in UMAP
FeaturePlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nFeature_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with \nnumber of detected genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by number of feature count or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.nFeature_RNA") + 
           ggtitle( "Map of cells by filter status \non nFeature_RNA value")
  )
}

# ..........................................................................................................
## @knitr dimreduc_ggplot_covariables_tSNE
# ..........................................................................................................

# Plot the cells percentage of ribosomal genes in tSNE
FeaturePlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "tsne"), 
             features = "percent.rb") +
              scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
              ggtitle( "Map of cells with level of\n percentage of ribosomal genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by %of ribosomal genes or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "tsne"),
                                  group.by = "outlier.percent.rb") + 
                 ggtitle( "Map of cells by filter status\n on %ribosomal gene")
      )
}

# Plot the cells percentage of mitochondrial genes in tSNE  
FeaturePlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "tsne"),
             features = "percent.mt") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with level of percentage of\n mitochondrial genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by %of mitochondrial genes or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "tsne"),
                  group.by = "outlier.percent.mt") + 
           ggtitle( "Map of cells by filter status\n on %mitochondrial gene")
  )
}

# Plot the cells RNA counts in tSNE
FeaturePlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "tsne"),
             features = "nCount_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with RNA counts")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by number of UMI count or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "tsne"),
                  group.by = "outlier.nCount_RNA") + 
           ggtitle( "Map of cells by filter status \non nCount_RNA value")
  )
}

# Plot the cells genes counts in tSNE
FeaturePlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "tsne"),
             features = "nFeature_RNA") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with \nnumber of detected genes")

# if QC EXPLORATION MODE is active, plot if the cells are excluded by number of feature count or not
if( QC_EXPLORATION_MODE == TRUE){
  print( DimPlot( seurat_obj, reduction = ifelse(exists("useReduction"), useReduction, "tsne"),
                  group.by = "outlier.nFeature_RNA") + 
           ggtitle( "Map of cells by filter status \non nFeature_RNA value")
  )
}