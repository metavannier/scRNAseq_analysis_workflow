## @knitr subset
# Load input
SEURAT_OBJECT = snakemake@input[["seurat_object"]]
#Load path
PATH <- getwd()
PATH <- gsub("/03_Script","", PATH)
# Load output

# Load parameters
MARKERGENE = snakemake@params[["markergene"]]
GENETOREMOVE = snakemake@params[["genetoremove"]]
FEATURES = snakemake@params[["features"]]
CLUSTER = snakemake@params[["cluster"]]
WT = snakemake@params[["wt"]]
TEST = snakemake@params[["test"]]
THRESHOLD = snakemake@params[["threshold"]]
FCCUTOFF  = snakemake@params[["FCcutoff"]]
PCCUTOFF  = snakemake@params[["pCutoff"]]

# Load seurat object generated after clustering
obj_seurat <- readRDS(SEURAT_OBJECT)
clusterlist = c(strsplit(CLUSTER, ",")[[1]])
featureslist = c(strsplit(FEATURES, ",")[[1]])

# Subset all idents by gene expression level
for (i in 1:length(clusterlist)) {
    markerpositive <- subset(x = obj_seurat, subset = RIPOR2 > 0, idents = clusterlist[i])
    markernegative <- subset(x = obj_seurat, subset = RIPOR2 == 0, idents = clusterlist[i])

    new.cluster.ids <- "RIPOR2+"
    names(new.cluster.ids) <- levels(markerpositive)
    markerpositive <- RenameIdents(markerpositive, new.cluster.ids)
    markerpositive <- AddMetaData(
        object = markerpositive,
        metadata = new.cluster.ids,
        col.name = "expression")

    new.cluster.ids <- "RIPOR2-"
    names(new.cluster.ids) <- levels(markernegative)
    markernegative <- RenameIdents(markernegative, new.cluster.ids)
    markernegative <- AddMetaData(
        object = markernegative,
        metadata = new.cluster.ids,
        col.name = "expression")

    obj_seurat_cluster <- merge(x = markerpositive, y = markernegative)
    ### If you need to remove genes for the DE analyses 
    counts <- GetAssayData(obj_seurat_cluster, assay = "RNA")
    counts <- counts[-(which(rownames(counts) %in% c(GENETOREMOVE))),]
    obj_seurat_cluster <- subset(obj_seurat_cluster, features = rownames(counts))
    de.markers <- FindMarkers(obj_seurat_cluster, ident.1 = "RIPOR2+", ident.2 = "RIPOR2-", logfc.threshold = THRESHOLD, test.use = TEST)
    de.file = paste(PATH,"/05_Output/04_diffexp_subset/differential_expression/",clusterlist[i],"_subcluster_DE.csv",sep = "")
    write.csv(x = de.markers, file = de.file, quote = FALSE)
    if (i > 1){
        obj_seurat_all_cluster <- merge(x = obj_seurat_all_cluster, y = obj_seurat_cluster)
    }else{
        obj_seurat_all_cluster <- obj_seurat_cluster
    }
}

## @knitr volcano

de_files <- dir("../05_Output/04_diffexp_subset/differential_expression/")
for (k in 1:length(de_files)){
    de_name <- gsub("_subcluster_DE.csv","", de_files[[k]])
    de_file = paste(PATH,"/05_Output/04_diffexp_subset/differential_expression/",de_files[[k]],sep = "")
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

    sdr <- subset(rlog_results, (rlog_results$avg_log2FC  < -FC) & (rlog_results$p_val_adj < p))
    write.table(as.data.frame(sdr), file= paste("../05_Output/04_diffexp_subset/de_significant/",de_name,"_signif-down-regulated.txt", sep=""), quote = FALSE, sep = "\t")

    sur <- subset(rlog_results, (rlog_results$avg_log2FC > FC) & (rlog_results$p_val_adj < p))
    write.table(as.data.frame(sur), file= paste("../05_Output/04_diffexp_subset/de_significant/",de_name,"_signif-up-regulated.txt", sep=""), quote = FALSE, sep = "\t")
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
    title = paste("Differential expression between ", MARKERGENE,"+ and ", MARKERGENE, "-",sep = "")
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
    pdf(paste(PATH,"/05_Output/04_diffexp_subset/volcanoplot/volcano_", de_name, ".pdf",sep = ""))
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

## @knitr violin

ViolinPlot <- VlnPlot(obj_seurat_all_cluster, features = featureslist, split.by = "expression", group.by = "cluster")
plot(ViolinPlot)

for (i in 1:length(featureslist)) {
    pdffile = paste(PATH,"/05_Output/04_diffexp_subset/violin_plot/",featureslist[i],"_violin_plot.pdf",sep = "")
    pdf(pdffile)
    vln = VlnPlot(obj_seurat_all_cluster, features = featureslist[i], split.by = "expression", group.by = "cluster")
    print(vln)
    dev.off()
}

## @knitr info
sessionInfo()
