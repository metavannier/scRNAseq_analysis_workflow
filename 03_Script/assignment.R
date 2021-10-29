## @knitr assignment

TSNE = snakemake@input[["tsne"]]
DEMUXLET = snakemake@input[["demuxlet"]]
TABDEMUXLET = snakemake@output[["tabdemuxlet"]]

# read in the barcodes
tsne <- fread(TSNE)
demuxlet <- fread(DEMUXLET)

# filter for the barcodes that we sampled
demuxlet<-subset(demuxlet, BARCODE %in% tsne$Barcode)
df <- data.frame(tsne1=tsne$"TSNE-1"[na.omit(match(demuxlet$BARCODE,tsne$Barcode))], tsne2=tsne$"TSNE-2"[na.omit(match(demuxlet$BARCODE,tsne$Barcode))], doublet=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[1]]}),
cell.type=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[2]]}))

# plot sample assignments 
#png(file="demuxlet_colored_tsne_geno_error_001_alpha_0_to_05.png",width=1200, height=1200, res=96)
ggplot(aes(tsne1,tsne2,color=cell.type),data=df)+geom_point()
#dev.off()

# plot doublet predictions
#png(file="demuxlet_doublets_tsne_geno_error_001_alpha_0_to_05.png",width=1200, height=1200, res=96)
ggplot(aes(tsne1,tsne2,color=doublet),data=df)+geom_point()
#dev.off() 

a<-subset(demuxlet, select=c("BARCODE","BEST"))
write.table(file=TABDEMUXLET, a ,sep="\t", quote=F, row.names=F, col.names=F)

#cat demuxlet_with_private_SNPs_geno_error_001_alpha_0_to_05.best.tsv | tr '-' '\t' | cut -f1,3- > Mix_MM_57_74_87_SOX10_24h_combined/outs/demuxlet_with_private_SNPs_geno_error_001_alpha_0_to_05.best.tsv