# #########################################
# This script reads and filters sc10x  data
# #########################################

# ..........................................................................................................
## @knitr loadData
# ..........................................................................................................

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

MIN_CELLS <- as.numeric(MIN_CELLS)
MIN_FEATURES <- as.numeric(MIN_FEATURES)
HTO <- unlist(strsplit(HTO, ","))

#........................................
# Read data
#........................................
seurat_data <- Read10X(paste0(in_data_dir,SAMPLE_ID,CELL_RANGER_COUNT_PATH))

# Need a condition if data are depultiplexed with cell ranger multi
if( RUN_DEMULTIPLEX == TRUE){
        seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = MIN_CELLS, min.features = MIN_FEATURES, project = SAMPLE_ID)
} else if (as.logical(MULTIPLEX) == TRUE) {
        seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`, min.cells = MIN_CELLS, min.features = MIN_FEATURES, project = SAMPLE_ID) 
        hto <- seurat_data$`Antibody Capture`
        hto <- hto[c(HTO), ]
        seurat_obj[["HTO"]] <- CreateAssayObject(counts = hto)
        seurat_obj <- subset(x = seurat_obj, subset = nCount_HTO > 0)
        seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")       
} else{
        seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`, min.cells = MIN_CELLS, min.features = MIN_FEATURES, project = SAMPLE_ID) 
}

seurat_obj$SampleID <- SAMPLE_ID

# ..........................................................................................................
## @knitr filterData_selection
# ..........................................................................................................

#........................................
#  Filter data
#........................................
FILTER_UMI_MIN <- as.numeric(FILTER_UMI_MIN)
FILTER_UMI_MAX <- as.numeric(FILTER_UMI_MAX)
FILTER_FEATURE_MIN <- as.numeric(FILTER_FEATURE_MIN)
FILTER_FEATURE_MAX <- as.numeric(FILTER_FEATURE_MAX)
FILTER_PERCENT_MT <- as.numeric(FILTER_PERCENT_MT)
FILTER_PERCENT_MT_MIN <- as.numeric(FILTER_PERCENT_MT_MIN)
FILTER_PERCENT_RB <- as.numeric(FILTER_PERCENT_RB)


if( QC_EXPLORATION_MODE == TRUE){
  cat("<p style='color:red;'><b>WARNING: QC EXPLORATION MODE is active : Cells will not be filtered but only marked as different</b></p>")
}

#........................................
# Mitochondrial and Ribosomal
#........................................
# Identify mitocondrial genes in matrix
mt.genes = grep( pattern = PATTERN_MT, x = rownames(GetAssayData(object = seurat_obj, slot = "counts")), value = TRUE, ignore.case = TRUE)
if(length(mt.genes)==0) 
{
warning( "No mitochondrial genes could be identified in this dataset.");
} else 
{
# Compute percentage of mitochondrial transcripts in each cell
percent.mt <- PercentageFeatureSet(seurat_obj, features = mt.genes)
# Add the mitocondrial gene percentage as meta information in the Seurat object
seurat_obj[["percent.mt"]] <- percent.mt
}

# Identify ribosomal genes in matrix
rb.genes = grep( pattern = PATTERN_RB, x = rownames(GetAssayData(object = seurat_obj, slot = "counts")), value = TRUE, ignore.case = TRUE)
if(length(rb.genes)==0) 
{
warning( "No ribosomal genes could be identified in this dataset.");
} else 
{
# Compute percentage of ribosomal transcripts in each cell
percent.rb <- PercentageFeatureSet(seurat_obj, features = rb.genes)
# Add the ribosomal gene percentage as meta information in the Seurat object
seurat_obj[["percent.rb"]] <- percent.rb
}

#...................................................................
# Identify cells that will be rejected based on specified thresholds
#...................................................................
# Reject cells based on UMI numbers
nUMI.drop = logical( length(Cells(seurat_obj)));
if( ! is.null( FILTER_UMI_MIN))
{
nUMI.drop = nUMI.drop | (seurat_obj[["nCount_RNA", drop=TRUE]] < FILTER_UMI_MIN);
}
if( ! is.null( FILTER_UMI_MAX))
{
nUMI.drop = nUMI.drop | (seurat_obj[["nCount_RNA", drop=TRUE]] > FILTER_UMI_MAX);
}
seurat_obj = AddMetaData( seurat_obj, metadata = nUMI.drop[ Cells( seurat_obj)], col.name = "outlier.nCount_RNA")

# Reject cells based on number of expressed genes
nGene.drop = logical( length(Cells(seurat_obj)));
if( ! is.null( FILTER_FEATURE_MIN))
{
nGene.drop = nGene.drop | (seurat_obj[["nFeature_RNA", drop=TRUE]] < FILTER_FEATURE_MIN);
}
if( ! is.null( FILTER_FEATURE_MAX))
{
nGene.drop = nGene.drop | (seurat_obj[["nFeature_RNA", drop=TRUE]] > FILTER_FEATURE_MAX);
}
seurat_obj = AddMetaData( seurat_obj, metadata = nGene.drop[ Cells( seurat_obj)], col.name = "outlier.nFeature_RNA")

# Identify cells with high and low percentage of mitocondrial genes
mito.drop = logical( length(Cells(seurat_obj)));
if( length(mt.genes) && (! is.null(FILTER_PERCENT_MT)))
{
mito.drop = (seurat_obj[["percent.mt", drop=TRUE]] > FILTER_PERCENT_MT);
}
if( ! is.null( FILTER_PERCENT_MT_MIN))
{
mito.drop = mito.drop | (seurat_obj[["percent.mt", drop=TRUE]] < FILTER_PERCENT_MT_MIN);
}
seurat_obj = AddMetaData( seurat_obj, metadata = mito.drop[ Cells( seurat_obj)], col.name = "outlier.percent.mt")


# Identify cells with low percentage of ribosomal genes
ribo.drop = logical( length(Cells(seurat_obj)));
if( length(rb.genes) && (! is.null(FILTER_PERCENT_RB)))
{
ribo.drop = (seurat_obj[["percent.rb", drop=TRUE]] < FILTER_PERCENT_RB);
}
seurat_obj = AddMetaData( seurat_obj, metadata = ribo.drop[ Cells( seurat_obj)], col.name = "outlier.percent.rb")

#..................................................................
# Plot distributions of #UMIs, #Genes, %Mito, and %Ribo among cells
#..................................................................
if(length( Cells(seurat_obj)) < PLOT_RASTER_NBCELLS_THRESHOLD)
{
        # Do interactive plot using plotly (can cause overhead when viewing many point)

        # Gather data to be visualized together (cell name + numID + metrics)
        cellsData = cbind( "Cell" = colnames(seurat_obj), # cell names from rownames conflict with search field in datatable
                        seurat_obj[[c( "SampleID", "nCount_RNA", "nFeature_RNA")]],
                        "percent.mt" = if(length( mt.genes)) as.numeric(format(seurat_obj[["percent.mt", drop = TRUE]], digits = 5)) else NULL,
                        "percent.rb" = if(length( rb.genes)) as.numeric(format(seurat_obj[["percent.rb", drop = TRUE]], digits = 5)) else NULL);

        # Create text to show under cursor for each cell
        hoverText = do.call(paste, c(Map( paste, 
                                        c( "", 
                                        "Cell ID: ", 
                                        "# UMIs: ", 
                                        "# Genes: ", 
                                        if(length( mt.genes)) "% Mito: ", 
                                        if(length( rb.genes)) "% Ribo: "), 
                                        cellsData, 
                                        sep=""), 
                        sep="\n"));

        panelWidth =  190; # 800 when using subplot
        panelHeight = 400;

        # Generate plotly violin/jitter panels for #umis, #genes, %mitochondrial, and %ribosomal stats
        lypanel_umis  = plotViolinJitter( cbind( cellsData, outliers = nUMI.drop), 
                                        xAxisFormula = ~as.numeric(1), 
                                        yAxisFormula = ~nCount_RNA, 
                                        pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                        hoverText = hoverText, 
                                        xTicklabels = "# UMIs", 
                                        thresholdHigh = FILTER_UMI_MAX, 
                                        thresholdLow = FILTER_UMI_MIN, 
                                        panelWidth = panelWidth, 
                                        panelHeight = panelHeight);

        lypanel_genes = plotViolinJitter( cbind( cellsData, outliers = nGene.drop), 
                                        xAxisFormula = ~as.numeric(1), 
                                        yAxisFormula = ~nFeature_RNA, 
                                        pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                        hoverText = hoverText, 
                                        xTicklabels = "# Genes", 
                                        thresholdHigh = FILTER_FEATURE_MAX, 
                                        thresholdLow = FILTER_FEATURE_MIN, 
                                        panelWidth = panelWidth, 
                                        panelHeight = panelHeight);

        lypanel_mitos = if(length(mt.genes)) plotViolinJitter( cbind( cellsData, outliers = mito.drop), 
                                                                xAxisFormula = ~as.numeric(1), 
                                                                yAxisFormula = ~percent.mt, 
                                                                pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                                                hoverText = hoverText, 
                                                                xTicklabels = "% Mito", 
                                                                thresholdHigh = FILTER_PERCENT_MT, 
                                                                thresholdLow = FILTER_PERCENT_MT_MIN, 
                                                                panelWidth = panelWidth, 
                                                                panelHeight = panelHeight) else NULL;


        lypanel_ribos = if(length(rb.genes)) plotViolinJitter( cbind( cellsData, outliers = ribo.drop), 
                                                                xAxisFormula = ~as.numeric(1), 
                                                                yAxisFormula = ~percent.rb, 
                                                                pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                                                hoverText = hoverText, 
                                                                xTicklabels = "% Ribo", 
                                                                thresholdLow = FILTER_PERCENT_RB, 
                                                                panelWidth = panelWidth, 
                                                                panelHeight = panelHeight) else NULL;

        # Set panels as a list and define plotly config
        panelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
        panelsList = lapply(panelsList, config, displaylogo = FALSE, 
                        toImageButtonOptions = list(format='svg'), 
                        modeBarButtons = list(list('toImage'), 
                                                list('zoom2d', 'pan2d', 'resetScale2d')));

        # Control layout using flex because subplot is limited regarding plot sizing and alignment
        # 'browsable' required in console, not in script/document
        browsable(
        div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
                div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
                lapply(panelsList, div, style = paste("flex : 0 0 auto; margin: 5px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
        )

} else
{
        # Do the plots as simple images with ggplot when having a lot of points

        # Generate plotly violin/jitter panels for #umis, #genes, %mitochondrial, and %ribosomal stats
        ggpanel_umis  = ggplot( cbind(seurat_obj[["nCount_RNA"]], outliers = nUMI.drop, seurat_obj[["SampleID"]]), aes( y = nCount_RNA)) +  
                        geom_violin(aes(x = 0.5), fill = "darkblue", col = "#AAAAAA", width = 0.9,alpha = 0.3) +
                        geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
                        geom_hline( yintercept = c(FILTER_UMI_MIN,FILTER_UMI_MAX), 
                                        color = c( "blue", "red")[!sapply( list( FILTER_UMI_MIN,FILTER_UMI_MAX), is.null)], 
                                        alpha = 0.5,
                                        size = 1) +
                        labs(x = "UMIs", y = "") + 
                        theme( axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank(),
                                axis.title.y = element_blank(),
                                legend.position = "none") +
                        scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));

        ggpanel_genes = ggplot( cbind(seurat_obj[["nFeature_RNA"]], outliers = nGene.drop, seurat_obj[["SampleID"]]), aes( y = nFeature_RNA)) +  
                        geom_violin(aes(x = 0.5), fill = "darkblue", col = "#AAAAAA", width = 0.9,alpha = 0.3) +
                        geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
                        geom_hline( yintercept = c(FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), 
                                        color = c( "blue", "red")[!sapply( list( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), is.null)], 
                                        alpha = 0.5,
                                        size = 1) +
                        labs(x = "Genes", y = "") +
                        theme( axis.text.x = element_blank(), 
                                axis.ticks.x = element_blank(),
                                axis.title.y = element_blank(),
                                legend.position = "none") +
                        scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));

        ggpanel_mitos = if(length(mt.genes)) ggplot( cbind(seurat_obj[["percent.mt"]], outliers = mito.drop, seurat_obj[["SampleID"]]), aes( y = percent.mt)) + 
                                                geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
                                                geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
                                                geom_hline( yintercept = c(FILTER_PERCENT_MT_MIN, FILTER_PERCENT_MT), 
                                                        color = c( "blue", "red")[!sapply( list( FILTER_PERCENT_MT_MIN, FILTER_PERCENT_MT), is.null)], 
                                                        alpha = 0.5,
                                                        size = 1) +
                                                labs(x = "% Mito", y = "") +
                                                theme( axis.text.x = element_blank(), 
                                                        axis.ticks.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        legend.position = "none") +
                                                scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444")) else NULL;

        ggpanel_ribos = if(length(rb.genes)) ggplot( cbind(seurat_obj[["percent.rb"]], outliers = ribo.drop, seurat_obj[["SampleID"]]), aes( y = percent.rb)) + 
                                                geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
                                                geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
                                                geom_hline( yintercept = FILTER_PERCENT_RB, 
                                                        color = "blue", 
                                                        alpha = 0.5,
                                                        size = 1) +
                                                labs(x = "% Ribo", y = "") +
                                                theme( axis.text.x = element_blank(), 
                                                        axis.ticks.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        legend.position = "none") +
                                                scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444")) else NULL;

        # Use patchwork library to assemble panels
        print( ggpanel_umis + ggpanel_genes + ggpanel_mitos + ggpanel_ribos + plot_layout( nrow = 1))
}

cat( "<br>Number of cells removed based on number of UMIs:", sum( nUMI.drop));
cat( "<br>Number of cells removed based on number of genes:", sum( nGene.drop));
if(exists( "mito.drop")) cat( "<br>Number of cells removed based on high percentage of mitochondrial transcripts:", sum( mito.drop));
if(exists( "ribo.drop")) cat( "<br>Number of cells removed based on low percentage of ribosomal transcripts:", sum( ribo.drop));
cat( "\n<br>\n");

# Identify cells to exclude as union of cells with low nb UMI, low nb expressed genes, high pct mito genes, low pct ribo genes
seurat_obj[["outlier"]] = nUMI.drop  | 
                     nGene.drop | 
                     ( if(exists( "mito.drop")) mito.drop else FALSE ) | 
                     ( if(exists( "ribo.drop")) ribo.drop else FALSE );

cat("<br><br>Removed cells after filters:", sum( unlist(seurat_obj[["outlier"]] )));
cat("<br>Remaining cells after filters:", sum( ! unlist(seurat_obj[["outlier"]] )));
cat("\n<br>\n");

#................................
# Record which cells got rejected
#................................
# Export the excluded cells to file
# write.table( data.frame( cellid = names(which(nUMI.drop))), 
#              file= file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_excluded_cells_nbUMI.txt")), 
#              quote = FALSE, 
#              row.names = FALSE, 
#              col.names = TRUE, 
#              sep="\t");

# write.table( data.frame( cellid = names(which(nGene.drop))), 
#              file= file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_excluded_cells_nbGene.txt")), 
#              quote = FALSE, 
#              row.names = FALSE, 
#              col.names = TRUE, 
#              sep="\t");

# if(exists( "mito.drop"))
# {
#   write.table( data.frame( cellid = names(which(mito.drop))), 
#                file= file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_excluded_cells_pctMito.txt")), 
#                quote = FALSE, 
#                row.names = FALSE, 
#                col.names = TRUE, 
#                sep="\t");
# }

# if(exists( "ribo.drop"))
# {
#   write.table( data.frame( cellid = names(which(ribo.drop))), 
#                file= file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_excluded_cells_pctRibo.txt")), 
#                quote = FALSE, 
#                row.names = FALSE, 
#                col.names = TRUE, 
#                sep="\t");
# }

# ..........................................................................................................
## @knitr filterData_summaryPlot
# ..........................................................................................................

#...................................................
# Plot dispersion of excluded and non-excluded cells
#...................................................
# number of genes and number of UMIs
ggplot( seurat_obj[[]][order( seurat_obj[["outlier"]]),], # Plot FALSE first and TRUE after
        aes( x = nFeature_RNA, 
             y = nCount_RNA, 
             color = outlier)) + 
  geom_point( size = 0.5) +
  geom_vline( xintercept = c( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), 
              linetype = 2, 
              color = c( "blue", "red")[!sapply( list( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), is.null)], 
              alpha = 0.5) +
  geom_hline( yintercept = c( FILTER_UMI_MIN, FILTER_UMI_MAX), 
              linetype = 2, 
              color = c( "blue", "red")[!sapply( list( FILTER_UMI_MIN, FILTER_UMI_MAX), is.null)], 
              alpha = 0.5) +
  scale_color_manual( values = c( "TRUE" = "#FF0000AA", "FALSE" = "#44444444")) + 
  labs( x = "# Genes", y = "# UMIs") +
  theme( legend.position = "none")

# Mitochondrial vs ribosomal distributions
if(exists( "mito.drop") && exists( "ribo.drop"))
{
  ggplot( seurat_obj[[]][order( seurat_obj[["outlier"]]),], # Plot FALSE first and TRUE after
          aes( x = percent.rb, 
               y = percent.mt, 
               color = outlier)) + 
    geom_point( size = 0.5) +
    geom_vline( xintercept = FILTER_PERCENT_RB, linetype = 2, color = "blue", alpha = 0.5) +
    geom_hline( yintercept = FILTER_PERCENT_MT, linetype = 2, color = "red", alpha = 0.5) +
    scale_color_manual( values = c( "TRUE" = "#FF0000AA", "FALSE" = "#44444444")) + 
    labs( x = "% Ribosomal genes", y = "% Mitochondrial genes") +
    theme( legend.position = "none")
}

# ..........................................................................................................
## @knitr filterData_filterObject
# ..........................................................................................................

# Filter the excluded cells in the Seurat object
seurat_obj_nonFiltered = seurat_obj;
if( QC_EXPLORATION_MODE == FALSE){
  seurat_obj = seurat_obj[ , ! seurat_obj[[ "outlier", drop=TRUE ]] ];
}

# Save the list of remaining cells after selection during loading, and #Genes, #UMIs, pctMito, pctRibo
# write.table( data.frame( cellid = Cells(seurat_obj)), 
#              file= file.path( OUTPUTDIR, STEP2, SAMPLE_ID, paste0( SAMPLE_ID, "_selected_cells.txt")), 
#              quote = FALSE, 
#              row.names = FALSE, 
#              col.names = TRUE, 
#              sep="\t");

# ..........................................................................................................
## @knitr normalizeData
# ..........................................................................................................

# Normalize the count data present in a given assay
seurat_obj = NormalizeData( object = seurat_obj,
                       normalization.method = NORM_METHOD,
                       scale.factor = NORM_SCALE_FACTOR,
                       verbose = .VERBOSE)

# # Scales and centers features in the dataset
# seurat_obj = ScaleData( object = seurat_obj,
#                    verbose = .VERBOSE)
