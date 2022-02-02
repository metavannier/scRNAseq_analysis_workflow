## @knitr violin
# Load input
SEURAT_OBJECT = snakemake@input[["seurat_object"]]
#Load path
PATH <- getwd()
PATH <- gsub("/03_Script","", PATH)
# Load output
VIOLINPLOT = snakemake@output[["violinplot"]]
UMAPFEATURE = snakemake@output[["umapfeature"]]
TSNEFEATURE = snakemake@output[["tsnefeature"]]
RIDGEFEATURE = snakemake@output[["ridgefeature"]]
HEATMAPFEATURE = snakemake@output[["heatmapfeature"]]
# Load parameters
FEATURES = snakemake@params[["features"]]
CELLMARKER = snakemake@params[["cellmarker"]]
WT = snakemake@params[["wt"]]
TEST = snakemake@params[["test"]]
THRESHOLD = snakemake@params[["threshold"]]
FCCUTOFF  = snakemake@params[["FCcutoff"]]
PCCUTOFF  = snakemake@params[["pCutoff"]]

# Load seurat object generated after clustering
obj_seurat <- readRDS(SEURAT_OBJECT)

featureslist = c(strsplit(FEATURES, ",")[[1]])

ViolinPlot <- VlnPlot(obj_seurat, features = featureslist)
plot(ViolinPlot)

for (i in 1:length(featureslist)) {
    pdffile = paste(PATH,"/05_Output/03_diffexp/violin_plot/",featureslist[i],"_violin_plot.pdf",sep = "")
    pdf(pdffile)
    vln = VlnPlot(obj_seurat, features = featureslist[i])
    print(vln)
    dev.off()
}

## @knitr umapfeature

umapfeature <- FeaturePlot(obj_seurat, features = featureslist, reduction = "umap")
plot(umapfeature)

for (i in 1:length(featureslist)) {
    pdffile = paste(PATH,"/05_Output/03_diffexp/umap_plot/",featureslist[i],"_umapfeature_plot.pdf",sep = "")
    pdf(pdffile)
    umapfeature <- FeaturePlot(obj_seurat, features = featureslist[i], reduction = "umap")
    print(umapfeature)
    dev.off()
}

## @knitr tsnefeature

tsnefeature <- FeaturePlot(obj_seurat, features = featureslist, reduction = "tsne")
plot(tsnefeature)

for (i in 1:length(featureslist)) {
    pdffile = paste(PATH,"/05_Output/03_diffexp/tsne_plot/",featureslist[i],"_tsnefeature_plot.pdf",sep = "")
    pdf(pdffile)
    tsnefeature <- FeaturePlot(obj_seurat, features = featureslist[i], reduction = "tsne")
    print(tsnefeature)
    dev.off()
}

## @knitr ridgefeature

ridgefeature <- RidgePlot(obj_seurat, features = featureslist)
plot(ridgefeature)

for (i in 1:length(featureslist)) {
    pdffile = paste(PATH,"/05_Output/03_diffexp/ridge_plot/",featureslist[i],"_ridgefeature_plot.pdf",sep = "")
    pdf(pdffile)
    ridgefeature <- RidgePlot(obj_seurat, features = featureslist[i])
    print(ridgefeature)
    dev.off()
}

## @knitr heatmapfeature

heatmapfeature <- DoHeatmap(object = obj_seurat, features = featureslist)
plot(heatmapfeature)
pdf(HEATMAPFEATURE)
print(heatmapfeature)
dev.off()

## @knitr testdiffexp

# Find differentially expressed features between WT and all the cellular type of interest together
de.markers <- FindMarkers(obj_seurat, ident.1 = WT, ident.2 = NULL, logfc.threshold = THRESHOLD, test.use = TEST)
de.file = paste(PATH,"/05_Output/03_diffexp/differential_expression_features/",WT,"_vs_AllCellMarker_DE.csv",sep = "")
write.csv(x = de.markers, file = de.file, quote = FALSE)

# Find differentially expressed features between WT and cellular type of interest
cellmarkerlist = c(strsplit(CELLMARKER, ",")[[1]])
for (i in 1:length(cellmarkerlist)) {
    de.file = paste(PATH,"/05_Output/03_diffexp/differential_expression_features/",WT,"_vs_",cellmarkerlist[i],"_DE.csv",sep = "")
    de.markers <- FindMarkers(obj_seurat, ident.1 = WT, ident.2 = cellmarkerlist[i], logfc.threshold = THRESHOLD, test.use = TEST)
    head(de.markers)
    write.csv(x = de.markers, file = de.file, quote = FALSE)
}

## @knitr volcano

de_files <- dir("../05_Output/03_diffexp/differential_expression_features/")
for (k in 1:length(de_files)){
    de_name <- gsub("_DE.csv","", de_files[[k]])
    de_file = paste(PATH,"/05_Output/03_diffexp/differential_expression_features/",de_files[[k]],sep = "")
    rlog_results <- read.table(de_file, sep = ",", dec = ".", header = TRUE, quote = "", skip = "")
    # Assign the first column as rownames
    rownames(rlog_results) <- rlog_results[,1]
    rlog_results[,1] <- NULL
    ## Volcano plot
    FC <- log2(FCCUTOFF)
    p <- PCCUTOFF
    # Red point in front of the other to visualize marker genes
    if(length(featureslist)!=0){
        for (i in 1:length(featureslist)) {
                gene_name=featureslist[[i]]
            if ((gene_name %in% row.names(rlog_results))){
                line <- rlog_results[match(gene_name, table = rownames(rlog_results)), ]
                rlog_results <- rlog_results[-match(gene_name, table = rownames(rlog_results)), ]
                rlog_results <- rbind(rlog_results,line)
            }
        }
    }

    keyvals <- rep('grey75', nrow(rlog_results))
    names(keyvals) <- rep('Non marker genes', nrow(rlog_results))
    # keyvals.shape <- rep(1, nrow(rlog_results))
    # names(keyvals.shape) <- rep('Non marker genes', nrow(rlog_results))

    # keyvals[which(abs(rlog_results$avg_log2FC) > FC & rlog_results$p_val_adj > p)] <- 'grey50'
    # names(keyvals)[which(abs(rlog_results$avg_log2FC) > FC & rlog_results$p_val_adj > p)] <- 'log2FoldChange'

    # keyvals[which(abs(rlog_results$avg_log2FC) < FC & rlog_results$p_val_adj < p)] <- 'grey25'
    # names(keyvals)[which(abs(rlog_results$avg_log2FC)  < FC & rlog_results$p_val_adj < p)] <- '-Log10Q'

    # keyvals[which(rlog_results$avg_log2FC < -FC & rlog_results$p_val_adj < p)] <- 'blue2'
    # names(keyvals)[which(rlog_results$avg_log2FC  < -FC & rlog_results$p_val_adj < p)] <- 'Signif. down-regulated'
    sdr <- subset(rlog_results, (rlog_results$avg_log2FC  < -FC) & (rlog_results$p_val_adj < p))
    write.table(as.data.frame(sdr), file= paste("../05_Output/03_diffexp/de_significant/",de_name,"_signif-down-regulated.txt", sep=""), quote = FALSE, sep = "\t")

    # keyvals[which(rlog_results$avg_log2FC > FC & rlog_results$p_val_adj < p)] <- 'red2'
    # names(keyvals)[which(rlog_results$avg_log2FC > FC & rlog_results$p_val_adj < p)] <- 'Signif. up-regulated'
    sur <- subset(rlog_results, (rlog_results$avg_log2FC > FC) & (rlog_results$p_val_adj < p))
    write.table(as.data.frame(sur), file= paste("../05_Output/03_diffexp/de_significant/",de_name,"_signif-up-regulated.txt", sep=""), quote = FALSE, sep = "\t")
    # Inlight the marker genes
    if(length(featureslist)!=0){
        for (i in 1:length(featureslist)) {
            gene_name=featureslist[[i]]
            if ((gene_name %in% row.names(rlog_results))){
                keyvals[which(row.names(rlog_results) == gene_name)] <- 'red2'
                names(keyvals)[which(row.names(rlog_results) == gene_name)] <- 'MarkerGenes'
            }
        }
    }
    unique(keyvals)
    unique(names(keyvals))
    # Title for the volcanoplot
    title = paste("Differential expression between ", de_name,sep = "")
    ## Volcanoplot with the adjusted p-value
    volcanoplot_padj <- EnhancedVolcano(rlog_results,
    lab = rownames(rlog_results),
    # selectLab = c('SOX10','SOX9','RIPOR2','MITF','DCT','MYC','FLT1','TNFRSF11B','XIRP2','ST3GAL1'),
    selectLab = featureslist,
    x = 'avg_log2FC',
    y = "p_val_adj",
    xlim = c(-5,5),
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    title = title,
    subtitle = "",
    axisLabSize = 12,
    titleLabSize = 15,
    pCutoff = p,
    FCcutoff = FC,
    pointSize = 1.5,
    labSize = 2,
    colCustom = keyvals,
    # shapeCustom = keyvals.shape,
    colAlpha = 1,
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 1.5,
    #DrawConnectors = TRUE,
    #widthConnectors = 0.2,
    colConnectors = 'grey50',
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = 'full',
    borderWidth = 0.7,
    borderColour = 'black')    
    plot(volcanoplot_padj)
    pdf(paste(PATH,"/05_Output/03_diffexp/volcanoplot/volcano_",de_name,".pdf",sep = ""))
    print(volcanoplot_padj)
    dev.off()
    # Tell if the list of interest genes have differential expression 
    if(length(featureslist)!=0){
        for (i in 1:length(featureslist)) {
            # adjusted p-value can be NA (>1) so condition doesn't work
            gene_name=featureslist[[i]]
            if ((gene_name %in% row.names(rlog_results))){
                foldc = rlog_results[which(row.names(rlog_results) == gene_name),"avg_log2FC"]
                pval = rlog_results[which(row.names(rlog_results) == gene_name),"p_val_adj"]

                if (!is.na(pval)){
                    if ((foldc < -FC) && (pval < p)) {
                        cat(paste("\n",gene_name," is differentially expressed (significantly down regulated).\n",sep=""))
                    }
                    if ((foldc > FC) && (pval < p)) {
                        cat(paste("\n",gene_name," is differentially expressed (significantly up regulated).\n",sep=""))
                    }
                }
            }
        }
    }
}

## @knitr info
sessionInfo()
