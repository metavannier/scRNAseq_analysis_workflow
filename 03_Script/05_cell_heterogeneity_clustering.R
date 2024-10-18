# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

# ..........................................................................................................
## @knitr heterogeneity_identifyClusters
# ..........................................................................................................

source(file.path(dirname(DIRECTORY), "03_Script/00_general_deps.R"))

#........................................
# Identify clusters
#........................................

# Computes the k.param nearest neighbors for a given dataset
seurat_obj <- FindNeighbors(object    = seurat_obj,
                            reduction = "pca",
                            dims      = 1:pcs,
                            verbose   = .VERBOSE);

# Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm
seurat_obj <- FindClusters(object             = seurat_obj,
                           resolution         = FINDCLUSTERS_RESOLUTION,
                           algorithm          = FINDCLUSTERS_ALGORITHM,
                           verbose            = .VERBOSE);

Idents(seurat_obj) = "seurat_clusters"

# Show number of cells in each cluster
clustersCount = as.data.frame( table( Cluster = seurat_obj[["seurat_clusters"]]), responseName = "CellCount");

# Define a set of colors for clusters (based on ggplot default)
clustersColor = hue_pal()( nlevels( Idents(seurat_obj)));
names( clustersColor) = levels( Idents(seurat_obj));

# Save cells cluster identity as determined with 'FindClusters'
# write.table( data.frame(seurat_obj[["SampleID"]], identity = Idents(seurat_obj)), 
#              file = file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_cellsClusterIdentity.tsv")), 
#              quote = FALSE, 
#              row.names = TRUE, 
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");

# # Also save cluster color attribution for reference
# # Save cells cluster identity as determined with 'FindClusters'
# write.table( clustersColor, 
#              file = file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_clustersColor.tsv")), 
#              quote = FALSE, 
#              row.names = TRUE, 
#              col.names = FALSE,
#              sep="\t");

# Create datatable
datatable( clustersCount,
           class = "compact",
           rownames = FALSE,
           colnames = c("Cluster", "Nb. Cells"),
           options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          columnDefs = list( # Center all columns
                                            list( targets = 0:(ncol(clustersCount)-1),
                                            className = 'dt-center')),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          paging = FALSE, # Disable pagination (show all)
                          processing = TRUE,
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE)) %>%

 # Add color from cluster
formatStyle( columns = "Cluster",
            color = styleEqual( names(clustersColor), clustersColor),
            fontWeight = 'bold')

# Compute a matrix of average expression value by cluster (for each gene)
geneExpByCluster = do.call( rbind, 
                            apply( as.matrix( GetAssayData(seurat_obj)), # expression values
                                   1,                                # by rows
                                   tapply,                           # apply by group
                                   INDEX = Idents(seurat_obj),       # clusters IDs
                                   mean,                             # summary function
                                   simplify = FALSE));               # do not auto coerce  

# Save it as 'tsv' file
# write.table( geneExpByCluster, 
#              file= file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_normExpressionByCluster.tsv")), 
#              quote = FALSE, 
#              row.names = TRUE, 
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");

# ..........................................................................................................
## @knitr heterogeneity_identifyClusters_splitStats
# ..........................................................................................................

#.....................................................................
# Show UMIs Genes Mitochondrial and Ribosomal content split by cluster
#.....................................................................
# Gather data to be visualized together (cell name + numID + metrics)
Idents(seurat_obj) = "seurat_clusters"
cellsData = cbind( "Cell" = colnames(seurat_obj), # cell names from rownames conflict with search field in datatable
                   seurat_obj[[c( "SampleID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mt" = if(length( mt.genes)) as.numeric(format(seurat_obj[["percent.mt", drop = TRUE]], digits = 5)) else NULL,
                   "percent.rb" = if(length( rb.genes)) as.numeric(format(seurat_obj[["percent.rb", drop = TRUE]], digits = 5)) else NULL,
                   "Batch" = seurat_obj[["orig.ident"]],
                   "Cluster" = Idents(seurat_obj));

# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste,
                                  c( "",
                                     "Cell ID: ",
                                     "# UMIs: ",
                                     "# Genes: ",
                                     if(length( mt.genes)) "% Mito: ",
                                     if(length( rb.genes)) "% Ribo: "),
                                  cellsData[-ncol( cellsData)],
                                  sep = ""),
                    sep = "\n"));

# Define size for panels (or assembled figure when using subplot)
panelWidth = 90 * nlevels(cellsData[["Cluster"]]);
panelHeight = 800;

# Generate plotly violin/jitter panels for #umis, #genes, and %mitochondrial stats
lypanel_umis  = plotViolinJitter( cellsData,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nCount_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# UMIs",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = panelWidth,
                                  panelHeight = panelHeight);

lypanel_genes = plotViolinJitter( cellsData,
                                  xAxisFormula = ~as.numeric( Cluster),
                                  yAxisFormula = ~nFeature_RNA,
                                  colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                  yAxisTitle = "# Genes",
                                  hoverText = hoverText,
                                  traceName = ~paste( "Cluster", Cluster),
                                  xTicklabels = levels(cellsData[["Cluster"]]),
                                  panelWidth = panelWidth,
                                  panelHeight = panelHeight);

lypanel_mitos = if(length( mt.genes)) plotViolinJitter( cellsData,
                                                          xAxisFormula = ~as.numeric( Cluster),
                                                          yAxisFormula = ~percent.mt,
                                                          colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                          yAxisTitle = "% Mito",
                                                          hoverText = hoverText,
                                                          traceName = ~paste( "Cluster", Cluster),
                                                          xTicklabels = levels(cellsData[["Cluster"]]),
                                                          panelWidth = panelWidth,
                                                          panelHeight = panelHeight) else NULL;

lypanel_ribos = if(length( rb.genes)) plotViolinJitter( cellsData,
                                                          xAxisFormula = ~as.numeric( Cluster),
                                                          yAxisFormula = ~percent.rb,
                                                          colorFormula = ~scales::alpha( clustersColor, 0.7)[as.character(Cluster)],
                                                          yAxisTitle = "% Ribo",
                                                          hoverText = hoverText,
                                                          traceName = ~paste( "Cluster", Cluster),
                                                          xTicklabels = levels(cellsData[["Cluster"]]),
                                                          panelWidth = panelWidth,
                                                          panelHeight = panelHeight) else NULL;

# Set panels as a list, define plotly config, and remove legend
panelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
panelsList = lapply( panelsList, config, displaylogo = FALSE,
                     toImageButtonOptions = list( format='svg'),
                     modeBarButtons = list( list('toImage'),
                                            list( 'zoom2d', 'pan2d', 'resetScale2d')));
                                        
# Group plotly violin/jitter panels so we can synchronise axes and use highlight on the full set
plotPanels = layout( subplot( panelsList,
                              nrows = 4,
                              shareX = TRUE,
                              titleY = TRUE),
                     xaxis = list(title = "Seurat Cluster",
                                  showgrid = TRUE,
                                  tickvals = seq(nlevels(cellsData[["Cluster"]])),
                                  ticktext = levels(cellsData[["Cluster"]])),
                     showlegend = FALSE, # Remove eventual legends (does not mix well with subplot and highlight)
                     autosize = TRUE);

# Control layout using flex
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
     div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
          div(plotPanels, style = paste("flex : 0 0 auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))

# ..........................................................................................................
## @knitr heterogeneity_dimReduc_with_clusters
# ..........................................................................................................

# Plot the map with clusters with ggplot
print(
DimPlot( seurat_obj, 
         reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
         group.by = "seurat_clusters",
         label = TRUE) + 
  ggtitle( "Map of cells with clusters")
)

# Plot non interactive version of the UMAP
if( exists( "CLUSTER_GROUP_LIST")){
  for( cluster_group_id in names( CLUSTER_GROUP_LIST)){
    
    # Get the cells in the cluster group
    Idents( seurat_obj) = "seurat_clusters"
    cells_in_group <- WhichCells( seurat_obj, idents = CLUSTER_GROUP_LIST[[ cluster_group_id]])
    
    # Highlight the cells in UMAP
    print(
      DimPlot( seurat_obj, 
               label=TRUE, label.size = 10,
               cells.highlight= list(cells_in_group), 
               cols.highlight = c( SAMPLE_COLOR[ SAMPLE_ALL]), 
               cols= "grey") +
           ggtitle( paste( 'UMAP of Cells in', cluster_group_id)) +
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none",
                  plot.margin = margin( 0, 0, 0, 0, "cm"),
                  plot.title = element_text( face = "bold",
                                             size = rel( 16/14),
                                             hjust = 0.5,
                                             vjust = 1,
                                             margin = margin( b = 7)))
    )
  }
}

# ..........................................................................................................
## @knitr heterogeneity_dimReduc_with_clusters_tSNE
# ..........................................................................................................

# Plot the map with clusters with ggplot
print(
DimPlot( seurat_obj, 
         reduction = ifelse(exists("useReduction"), useReduction, "tsne"), 
         group.by = "seurat_clusters",
         label = TRUE) + 
  ggtitle( "Map of cells with clusters")
)

# Plot non interactive version of the UMAP
if( exists( "CLUSTER_GROUP_LIST")){
  for( cluster_group_id in names( CLUSTER_GROUP_LIST)){
    
    # Get the cells in the cluster group
    Idents( seurat_obj) = "seurat_clusters"
    cells_in_group <- WhichCells( seurat_obj, idents = CLUSTER_GROUP_LIST[[ cluster_group_id]])
    
    # Highlight the cells in UMAP
    print(
      DimPlot( seurat_obj, 
               label=TRUE, label.size = 10,
               cells.highlight= list(cells_in_group), 
               cols.highlight = c( SAMPLE_COLOR[ SAMPLE_ALL]), 
               cols= "grey") +
           ggtitle( paste( 'UMAP of Cells in', cluster_group_id)) +
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none",
                  plot.margin = margin( 0, 0, 0, 0, "cm"),
                  plot.title = element_text( face = "bold",
                                             size = rel( 16/14),
                                             hjust = 0.5,
                                             vjust = 1,
                                             margin = margin( b = 7)))
    )
  }
}

# ..........................................................................................................
## @knitr save_RDS
# ..........................................................................................................

# Save each seurat object
saveRDS(seurat_obj, file = file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0(SAMPLE_ID,"_filtered_seurat_object.rds")))

# # save a tsv file for each sample
# write.table( seurat_obj[[]], 
#              file= file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_filtered_seurat_object.tsv")), 
#              quote = FALSE, 
#              row.names = TRUE, 
#              col.names = NA, # Add a blank column name for row names (CSV convention)
#              sep="\t");

# # Save the count matrix as csv file 
# count_matrix <- GetAssayData(object = seurat_obj[["RNA"]], slot = "counts")
# count_matrix <- as.matrix(count_matrix)
# count_matrix <- t(count_matrix)
# count_matrix <- as.data.frame(count_matrix)

# fwrite(x = count_matrix, file = file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_count_matrix.csv")), row.names = TRUE)

# # Save the normalize matrix as csv file
# normalized_matrix <- GetAssayData(object = seurat_obj[["RNA"]], slot = "data")
# normalized_matrix <- as.matrix(normalized_matrix)
# normalized_matrix <- t(normalized_matrix)
# normalized_matrix <- as.data.frame(normalized_matrix)

# fwrite(x = normalized_matrix , file = file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_normalized_matrix.csv")), row.names = TRUE)
