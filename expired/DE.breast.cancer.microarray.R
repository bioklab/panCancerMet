#REF: Transcriptional Profiling of Breast Cancer Metastases Identifies Liver Metastasis–Selective Genes Associated with Adverse Outcome in Luminal A Primary Breast Cancer

#GSE46141
#GSE14018
#GSE14020




GSE46141_Probeset_Annotation_20100204           <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/GSE46141_Probeset_Annotation_20100204.txt", stringsAsFactors=FALSE)
rownames(GSE46141_Probeset_Annotation_20100204) <- GSE46141_Probeset_Annotation_20100204$probesetid



require(GEOquery)
gse <- getGEO("GSE46141", GSEMatrix = FALSE)
expr.matrix <- foreach(gsm.obj = gse@gsms,.combine='cbind') %do% {
    gsm.obj@dataTable@table 
}
flag                  <- grepl(x=colnames(expr.matrix),pattern = 'VALUE')
rownames(expr.matrix) <- expr.matrix$ID_REF
expr.matrix           <- expr.matrix[,flag]
colnames(expr.matrix) <- names(gse@gsms)
expr.matrix           <- as.matrix(expr.matrix)


biospy.site <- foreach(gsm.obj = gse@gsms,.combine='c') %do% {
    gsm.obj@header$source_name_ch1
}
sample.annotation           <- data.frame(gsm.id = colnames(expr.matrix),biospy.site=biospy.site)
rownames(sample.annotation) <- sample.annotation$gsm.id



probe.set <- GSE46141_Probeset_Annotation_20100204$probesetid[GSE46141_Probeset_Annotation_20100204$GeneSymbol %in% 'ESR1']
pca.rs    <- prcomp(expr.matrix[probe.set,] %>% t)
plot(pca.rs$x[,1:2])
pc1       <- pca.rs$x[,1]
sample.annotation$subtype                        <- 'basal'
sample.annotation[names(pc1)[pc1 < 0],'subtype'] <- 'ER'

flag <- sample.annotation$subtype != 'basal' & sample.annotation$biospy.site %in% c('liver','breast')

s    <- sample.annotation[flag,]
eset <- expr.matrix[,flag]
s$biospy.site <- as.character(s$biospy.site)
design <- model.matrix(data = s, object= ~ biospy.site)

fit <- lmFit(eset, design)
ebayes <- eBayes(fit)
p.value.vec <- (ebayes[["p.value"]][,2]) %>% sort
tmp <- data.frame(gene.symbol=GSE46141_Probeset_Annotation_20100204[names(p.value.vec),'GeneSymbol'],p.value=p.value.vec)
tmp <- tmp[tmp$gene.symbol != '',]

combine.p.value.fisher <- function(df) {
  fisher.statistics     <- -2 * sum(log(df$p.value))
  p.value               <- pchisq(q=fisher.statistics,df = 2 * nrow(df),lower.tail = FALSE)
  rs                    <- data.frame(p.value = p.value)
}

cp.p.value.df <- ddply(tmp,.(gene.symbol),combine.p.value.fisher)
cp.p.value.df <- cp.p.value.df[order(cp.p.value.df$p.value),]
cp.p.value.df$fdr <- p.adjust(cp.p.value.df$p.value,method='fdr')
de.gene <- cp.p.value.df$gene.symbol[cp.p.value.df$fdr < 0.05] %>% as.character()

rownames(cp.p.value.df) <- cp.p.value.df$gene.symbol

df <- cp.p.value.df[basal.up.gene.df$SYMBOL %>% as.character(),]
df <- df[order(df$p.value),]
#fdr.vec <- p.adjust(p.value.vec,method='fdr')
#fdr.vec <- sort(fdr.vec)


de.res.liver.vs.breast.up.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = de.res.liver.vs.breast.up.gene)
de.res.liver.vs.breast.dn.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = de.res.liver.vs.breast.dn.gene)



GSE46141_Probeset_Annotation_20100204[names(fdr.vec)[fdr.vec < 0.1],] %>% View



# m.matrix                    <- expr.matrix[,sample.annotation$gsm.id[sample.annotation$biospy.site == 'liver'] %>% as.character()] %>% as.matrix()
# p.matrix                    <- expr.matrix[,sample.annotation$gsm.id[sample.annotation$biospy.site == 'breast'] %>% as.character()] %>% as.matrix()
# 
# 
# c.gene <- intersect(GSE46141_Probeset_Annotation_20100204$GeneSymbol,CCLE.rna.seq.marker.gene.1000.annotation$SYMBOL)
# probe.set <- GSE46141_Probeset_Annotation_20100204$probesetid[GSE46141_Probeset_Annotation_20100204$GeneSymbol %in% c.gene]
# dist.matrix <- 1- cor(cbind(m.matrix,p.matrix)[probe.set,] , method='spearman')
# tsne.rs <- Rtsne(as.dist(dist.matrix),perplexity = 4)
# plot(tsne.rs$Y)
# rownames(tsne.rs$Y) <- c(colnames(m.matrix),colnames(p.matrix))
# 
# 
# p.value.vec <- sapply(1:nrow(m.matrix),function(i) wilcox.test(m.matrix[i,],p.matrix[i,])$p.value   )
# names(p.value.vec) <- rownames(m.matrix)
# p.value.vec <- sort(p.value.vec)
# 
# 
# lumb.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv",  stringsAsFactors=FALSE)$x
# basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
# lumb.up.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = lumb.up.gene)
# basal.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = basal.up.gene)
# 
# probe.set <- GSE46141_Probeset_Annotation_20100204$probesetid[GSE46141_Probeset_Annotation_20100204$GeneSymbol %in% lumb.up.gene.df$SYMBOL]
# 
# 
# 
# pam50_annotation <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/pam50_annotation.txt", stringsAsFactors=FALSE)
# probe.set <- GSE46141_Probeset_Annotation_20100204$probesetid[GSE46141_Probeset_Annotation_20100204$GeneSymbol %in% pam50_annotation$GeneName]
# probe.set <- GSE46141_Probeset_Annotation_20100204$probesetid[GSE46141_Probeset_Annotation_20100204$GeneSymbol %in% 'ESR1']
# 
# pca.rs <- prcomp(expr.matrix[probe.set,] %>% t)
# plot(pca.rs$x[,1:2])
# pc1 <- pca.rs$x[,1]
# sample.annotation$subtype <- 'basal'
# sample.annotation[names(pc1)[pc1 < 0],'subtype'] <- 'ER'
# 
# m.matrix                    <- expr.matrix[,sample.annotation$gsm.id[sample.annotation$biospy.site == 'liver' & sample.annotation$subtype != 'basal'] %>% as.character()] %>% as.matrix()
# p.matrix                    <- expr.matrix[,sample.annotation$gsm.id[sample.annotation$biospy.site == 'breast'& sample.annotation$subtype != 'basal'] %>% as.character()] %>% as.matrix()
# p.value.vec <- sapply(1:nrow(m.matrix),function(i) wilcox.test(m.matrix[i,],p.matrix[i,])$p.value   )
# names(p.value.vec) <- rownames(m.matrix)
# p.value.vec <- sort(p.value.vec)
# 
# probe.set <- GSE46141_Probeset_Annotation_20100204$probesetid[GSE46141_Probeset_Annotation_20100204$GeneSymbol %in% basal.up.gene.df$SYMBOL]
# hist(p.value.vec[probe.set])

# require(affy)
# require(oligo)
# cel.files <- list.celfiles('/Users/liuke/Project/BreastCancerMetaPotenial/client-side/Data/GSE46141_RAW/',full.names=TRUE)
# rawData <- read.celfiles(cel.files)
# 
# require(makecdfenv) #https://www.biostars.org/p/67400/    very useful!
# pkgpath <- tempdir()
# make.cdf.package("GPL10379_HuRSTA-2a520709_custom_MMPM.cdf", cdf.path="/Users/liuke/Project/BreastCancerMetaPotenial/client-side/Data/GSE46141_RAW/", compress=FALSE, species = "Homo_sapiens", package.path = pkgpath)
# dir(pkgpath)
# 
# require(gpl10379hursta2a520709custommmpmcdf)
# 
# data        <- ReadAffy(filenames = cel.files,cdfname ='gpl10379hursta2a520709custommmpmcdf' )
# data.norm   <- rma(d,normalize = TRUE,background = TRUE)
# expr.matrix <- exprs(d.norm)


## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11078
require(GEOquery)
gse <- getGEO("GSE11078", GSEMatrix = FALSE)
expr.matrix <- foreach(gsm.obj = gse@gsms,.combine='cbind') %do% {
  gsm.obj@dataTable@table 
}
flag                  <- grepl(x=colnames(expr.matrix),pattern = 'VALUE')
rownames(expr.matrix) <- expr.matrix$ID_REF
expr.matrix           <- expr.matrix[,flag]
colnames(expr.matrix) <- names(gse@gsms)
expr.matrix           <- as.matrix(expr.matrix)

