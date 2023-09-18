# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

## @knitr ldr

# Scales and centers features in the dataset
seurat_obj = ScaleData( object = seurat_obj,
                   verbose = .VERBOSE)

# Perform linear dimensional reduction
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
VizDimLoadings(seurat_obj, dims = 1:DIMS, reduction = "pca")
# DimPlot(seurat_obj, reduction = "pca", group.by = "orig.ident")
# DimPlot(seurat_obj, reduction = "pca", group.by = "categorie")
DimHeatmap(seurat_obj, dims = 1:DIMS, cells = 500, balanced = TRUE)

## @knitr elbowplot
ElbowPlot(seurat_obj, ndims = DIMS)

# ..........................................................................................................
## @knitr heterogeneity_pca
# ..........................................................................................................

# Compute PCA on selected variable genes
nbPC = PCA_NPC
if(PCA_NPC > length(Cells(seurat_obj)))
{
  warning( paste0( "Number of cells in object (", length(Cells(seurat_obj)), ") smaller than requested number of PCs (", PCA_NPC,"), setting lower PC number..." ))
  nbPC = length(Cells(seurat_obj))
}           
seurat_obj <- RunPCA( object = seurat_obj,
                 features = VariableFeatures( seurat_obj),
                 npcs     = nbPC,
                 verbose  = .VERBOSE);

# Plot PCA, highlighting seurat clusters for combination of dimensions
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  print( DimPlot( object = seurat_obj, reduction = "pca", dims = dims, group.by = "orig.ident") +
           theme( legend.position = "none"));                              
}));

# ..........................................................................................................
## @knitr heterogeneity_pca_umisCounts
# ..........................................................................................................

# Same PCA plots but highlight UMI counts
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  print( FeaturePlot( seurat_obj, feature = "nCount_RNA", reduction = "pca", dims = dims) +
           theme( legend.position = "none",
                  plot.margin = margin(0, 0, 0, 0, "cm"),
                  plot.title = element_blank()));
}));

# ..........................................................................................................
## @knitr heterogeneity_pca_genesCounts
# ..........................................................................................................

# Same PCA plots but highlight feature counts
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  print( FeaturePlot( seurat_obj, feature = "nFeature_RNA", reduction = "pca", dims = dims) +
           theme( legend.position = "none",
                  plot.margin = margin(0, 0, 0, 0, "cm"),
                  plot.title = element_blank()));
}));

# ..........................................................................................................
## @knitr heterogeneity_pca_correlations
# ..........................................................................................................

# Isolate the PCA location of cells in first dimensions, together with UMI and genes counts for correlation analysis (makes use of cbind recycling to repeat values for each stacked PC)
relationToPC = suppressWarnings( cbind( stack( as.data.frame( Embeddings( seurat_obj, reduction = "pca")[ rownames( seurat_obj[[]]), paste0( "PC_", 1:PCA_PLOTS_NBDIMS) ])),
                                        seurat_obj[[ c( "nCount_RNA", "nFeature_RNA") ]],
                                        Cluster = Idents( seurat_obj)));

# Plot relationship of UMIs and genes counts with PCs (combine plots using '/' from 'patchwork' lib)
print( (ggplot( data = relationToPC, aes( x = values, y = nCount_RNA)) +
          facet_wrap( ~ind) +
          stat_binhex( bins = 60) +
          #geom_point( aes(col = Cluster), alpha = 0.5) +
          geom_smooth( method = 'lm') +
          stat_cor( method = "spearman") +
          ylab( "# UMIs") +
          theme( axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 plot.margin = unit( c( 1, 1, -0.5, 0.5), "lines")))
       /
         (ggplot( data = relationToPC, aes( x = values, y = nFeature_RNA)) +
            facet_wrap( ~ind) +
            stat_binhex( bins = 60) +
            #geom_point( aes(col = Cluster), alpha = 0.5) +
            geom_smooth( method = 'lm') +
            stat_cor( method = "spearman") +
            xlab( "PC values") +
            ylab( "# Genes")));

# ..........................................................................................................
## @knitr heterogeneity_pca_loadings
# ..........................................................................................................

# Plot PCA loadings
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  namesPC=paste0( "PC_", dims);
  # Get the loading values for concerned PCs
  loadingsMatrix = Loadings( seurat_obj, reduction = "pca")[ , namesPC ];
  # Sort features by average absolute value and convert as DF with features names as column
  loadingsMatrix = head( loadingsMatrix[ order( apply( loadingsMatrix, 1, function(x){ mean( abs( x)) }), decreasing = TRUE), ], PCA_PLOTS_NBFEATURES);
  loadingsDF = data.frame( loadingsMatrix, features = rownames( loadingsMatrix));
  
  # Define symmetric and consistent axes for group of plots
  axesLimit = max( abs( loadingsMatrix));
  
  # Plot arrows and features name
  print( ggplot( data = loadingsDF, aes_( x = as.name( namesPC[1]), y = as.name( namesPC[2]))) +
           coord_cartesian( xlim = c( -axesLimit, axesLimit), ylim = c( -axesLimit, axesLimit)) +
           geom_text_repel( aes( label = features), max.iter = 10000) + #, max.overlaps = Inf
           geom_segment( x = 0 , y = 0, aes_( xend = as.name(namesPC[1]), yend = as.name(namesPC[2])), col = "#00000044", arrow = arrow( length = unit( 2, "mm"))));
}));
