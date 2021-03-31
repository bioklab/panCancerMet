require(ggplot2)
require(CePa)
require(dplyr)
require(ComplexHeatmap)
require(RColorBrewer)
require(circlize)
require(segmented)
source('client-side/code/Manuscript/ggplot.style.R')
load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')


################################################################################   
#----------------------------- Figure 1 ---------------------------------------#
################################################################################  


################################################################################################
# Fig1a: cartoon to explain how the non-malignant cells influence DE analysis perfromed on bulk RNAseq data, made by illustrator (cartoon.pdf)
################################################################################################


################################################################################################
# Fig1b: scatter plot to show the relationship between log2FC values derived from differnt DE analysis comparisions
################################################################################################
c.gene        <- intersect(rownames(de.res.liver.vs.breast.lumb),rownames(de.res.metastasis.liver.vs.breast.lumb))
x             <- de.res.liver.vs.breast.lumb[c.gene,'log2FoldChange']
y             <- de.res.metastasis.liver.vs.breast.lumb[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=de.res.liver.vs.breast.lumb[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.lumb[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/lumb.log2FC.MP.and.normal.pdf',width = 20,height=20)


################################################################################################
# Fig 1c: an example of outlier DE CSF3R 
################################################################################################
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')

TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 

MET500.sample  <- c(MET500.breast.cancer.polyA.LumB.sample)
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]

TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample

df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
f1             <- df$condition == 'MET500'
f2             <- df$condition == 'TCGA'
g              <- 'ENSG00000119535' # An outlier example
df$expr        <- c(MET500.log2.fpkm.matrix[g,MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix[g,TCGA.sample])  
wilcox.test.p.value <- wilcox.test(df$expr[f1],df$expr[f2])$p.value

ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=2) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/outlier.example.pdf',width = 20,height=20)

de.res.metastasis.liver.vs.breast.lumb[g,'pvalue']
wilcox.test.p.value


################################################################################################
# Fig 1d: the pipeline figure made illustrator (pipeline.pdf)
################################################################################################








################################################################################   
#----------------------------- Figure S1 ---------------------------------------#
################################################################################   

################################################
# FigS1a: expression pattern of the highly-upregulated 75 genes across multiple tissues
################################################
median.tpm.matrix <- read.gct(file='client-side/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
get.gene.id <- function(x) {
  strsplit(x=x,split='\\.') %>% unlist %>% head(1)
}
rownames(median.tpm.matrix) <- sapply(rownames(median.tpm.matrix),get.gene.id)


up.gene   <- rownames(de.res.metastasis.liver.vs.breast.lumb)[de.res.metastasis.liver.vs.breast.lumb$log2FoldChange > 5  & de.res.metastasis.liver.vs.breast.lumb$padj < 0.05]
up.gene   <- intersect(up.gene,median.tpm.matrix %>% rownames)
dev.off()
pdf(file = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/lumb.up.gene.expression.pattern.heatmap.pdf',width = 20,height=20)
tmp           <- apply((median.tpm.matrix[up.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
Heatmap(tmp, col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = FALSE,show_heatmap_legend = FALSE)
dev.off()



draw.color.bar <- function(value.range, color.range){
  f <- colorRamp2(value.range,color.range)
  x <- sapply(value.range, f)
  plot(c(0, 40*length(x)), c(0, 100), type = "n")
  i <- 0:(length(x)-1)
  rect(i*40, 0, (i+1)*40, 100, col=x)
}
pdf(file = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/color.bar.pdf',width = 20,height=20)
draw.color.bar(value.range = seq(from=-5,to=5,by=0.5),
               color.range = colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev
)
dev.off()




# require(data.table)
# require(org.Hs.eg.db)
# require(AnnotationDbi)
# require(clusterProfiler)
# 
# up.gene   <- rownames(de.res.metastasis.liver.vs.breast.lumb)[de.res.metastasis.liver.vs.breast.lumb$log2FoldChange > 5  & de.res.metastasis.liver.vs.breast.lumb$padj < 0.05]
# df <- fread(input = 'client-side/Data/tissue.enriched.gene.csv') %>% as.data.frame
# liver.enriched.gene <- df$ensg_id[df$`enriched tissue` == 'liver']
# intersect(up.gene,liver.enriched.gene) %>% length
# up.gene %>% length
# (intersect(up.gene,liver.enriched.gene) %>% length) / (up.gene %>% length)
# 
# up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene) 
# GO.rs <- enrichGO(gene=up.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
# 
# compute.jaccard.index  <- function(x,y){
#   ( intersect(x,y) %>% length )/( unique(c(x,y))   %>% length )
#   
# }
# get.gene.list <- function(x) {
#   strsplit(x,split='\\/')[[1]] %>% unlist    
# }
# 
# 
# 
# enriched.term.filtering <- function(x) {
#   result    <- x %>% dplyr::filter(Count >=5 & qvalue <= 0.2 & pvalue<=0.001) %>% dplyr::arrange(pvalue) 
#   if(nrow(result) == 0){
#     return(NA)    
#   }
#   gene.str  <- result$geneID %>% unique
#   idx.vec   <- sapply(gene.str,function(x) which(result$geneID == x) %>% min  )
#   result    <- result[idx.vec,]
#   gene.list <- lapply(result$geneID,get.gene.list)
#   sim.matrix <- matrix(data=0,nrow=nrow(result),ncol=nrow(result))
#   for(i in 1:nrow(sim.matrix)){
#     for(j in 1:nrow(sim.matrix)){
#       sim.matrix[i,j] <- compute.jaccard.index(gene.list[[i]],gene.list[[j]])    
#     }
#   }
#   diag(sim.matrix) <- 0
#   
#   
#   bit.vec   <- rep(1,nrow(result))
#   for(i in 1:length(bit.vec)) {
#     if(bit.vec[i] == 1){
#       cover.index <- which(sim.matrix[i,] >= 0.4)
#       bit.vec[cover.index] <- 0
#     }        
#   }
#   result[bit.vec == 1,]
#   
# }
# GO.filtered.rs <- enriched.term.filtering(GO.rs)

################################################################################################
# FigS1b: show the confounding effects of liver cells in her2 subtype with MET500 dataset
################################################################################################
c.gene        <- intersect(rownames(de.res.liver.vs.breast.her2),rownames(de.res.metastasis.liver.vs.breast.her2))
x             <- de.res.liver.vs.breast.her2[c.gene,'log2FoldChange']
y             <- de.res.metastasis.liver.vs.breast.her2[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=de.res.liver.vs.breast.her2[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.her2[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/her2.log2FC.MP.and.normal.pdf',width = 20,height=20)


up.gene   <- rownames(de.res.metastasis.liver.vs.breast.her2)[de.res.metastasis.liver.vs.breast.her2$log2FoldChange > 5  & de.res.metastasis.liver.vs.breast.her2$padj < 0.05]
up.gene   <- intersect(up.gene,median.tpm.matrix %>% rownames)
dev.off()
pdf(file = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/her2.up.gene.expression.pattern.heatmap.pdf',width = 20,height=20)
tmp           <- apply((median.tpm.matrix[up.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
Heatmap(tmp, col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = FALSE,show_heatmap_legend = FALSE)
dev.off()

################################################################################################
# FigS1c : show the confounding effects of liver cells in basal subtype with MET500 dataset
################################################################################################
c.gene        <- intersect(rownames(de.res.liver.vs.breast.basal),rownames(de.res.metastasis.liver.vs.breast.basal))
x             <- de.res.liver.vs.breast.basal[c.gene,'log2FoldChange']
y             <- de.res.metastasis.liver.vs.breast.basal[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=de.res.liver.vs.breast.basal[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.basal[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/basal.log2FC.MP.and.normal.pdf',width = 20,height=20)



up.gene   <- rownames(de.res.metastasis.liver.vs.breast.basal)[de.res.metastasis.liver.vs.breast.basal$log2FoldChange > 5  & de.res.metastasis.liver.vs.breast.basal$padj < 0.05]
up.gene   <- intersect(up.gene,median.tpm.matrix %>% rownames)
dev.off()
pdf(file = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/basal.up.gene.expression.pattern.heatmap.pdf',width = 20,height=20)
tmp           <- apply((median.tpm.matrix[up.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
Heatmap(tmp, col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = FALSE,show_heatmap_legend = FALSE)
dev.off()


#############################################################
# FigS1d: use another dataset GSE58708 to show the same pattern
#############################################################
load('client-side/output/validation.of.confounding.R.output/validation.of.confounding.RData')

c.gene        <- intersect(rownames(SRP043470.de.res.liver.vs.breast),rownames(SRP043470.de.res.metastasis.liver.vs.breast))
x             <- SRP043470.de.res.liver.vs.breast[c.gene,'log2FoldChange']
y             <- SRP043470.de.res.metastasis.liver.vs.breast[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=SRP043470.de.res.liver.vs.breast[c.gene,'log2FoldChange'],y=SRP043470.de.res.metastasis.liver.vs.breast[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/SRP043470.log2FC.MP.and.normal.pdf',width = 20,height=20)



up.gene   <- rownames(SRP043470.de.res.liver.vs.breast)[SRP043470.de.res.liver.vs.breast$log2FoldChange > 5  & SRP043470.de.res.liver.vs.breast$padj < 0.05]
up.gene   <- intersect(up.gene,median.tpm.matrix %>% rownames)
dev.off()
pdf(file = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/SRP043470.up.gene.expression.pattern.heatmap.pdf',width = 20,height=20)
tmp           <- apply((median.tpm.matrix[up.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
Heatmap(tmp, col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = FALSE,show_heatmap_legend = FALSE)
dev.off()





################################################################################################
# Fig S1e: boxplot of CD8A expression
################################################################################################
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')

TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 


MET500.sample  <- c(MET500.breast.cancer.polyA.LumB.sample)
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]

TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample


df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
f1             <- df$condition == 'MET500'
f2             <- df$condition == 'TCGA'
g              <- 'ENSG00000153563' # CD8A
df$expr        <- c(MET500.log2.fpkm.matrix[g,MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix[g,TCGA.sample])  
wilcox.test.p.value <- wilcox.test(df$expr[f1],df$expr[f2])$p.value

ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=2) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/CD8A.example.pdf',width = 20,height=20)

de.res.metastasis.liver.vs.breast.lumb[g,'pvalue']
wilcox.test.p.value




##############################################################################
#---------------- Figure S2 ------------------------------------------------#
##############################################################################



##############################################################################
#Figure S2a:boxplot of MYC
##############################################################################

load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')


TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 


MET500.sample  <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample

MYC <- 'ENSG00000136997'

df1 <- data.frame(expr=MET500.log2.fpkm.matrix[MYC,MET500.sample],condition='MET500')
df2 <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[MYC,TCGA.sample],condition='TCGA')
df <- rbind(df1,df2)
df$condition <- factor(df$condition,levels = c('TCGA','MET500'))
ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=3) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') 

ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/MYC.DE.pdf',width = 20,height=20)




##############################################################################
#Figure S2b,S2c,S2d : explore the expression pattern of identified tumor-intrinsic DE genes
##############################################################################

lumb.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv",  stringsAsFactors=FALSE)$x
basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
her2.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv",  stringsAsFactors=FALSE)$x
lumb.dn.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.dn.csv",  stringsAsFactors=FALSE)$x
her2.dn.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.dn.csv",  stringsAsFactors=FALSE)$x
basal.dn.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.dn.csv", stringsAsFactors=FALSE)$x

lumb.gene  <- intersect(c(lumb.up.gene,lumb.dn.gene),  rownames(median.tpm.matrix))
her2.gene  <- intersect(c(her2.up.gene,her2.dn.gene),  rownames(median.tpm.matrix))
basal.gene <- intersect(c(basal.up.gene,basal.dn.gene),rownames(median.tpm.matrix))


tmp           <- apply((median.tpm.matrix[lumb.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
flag          <- which(rowSums(tmp) > 0)
pdf(file = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/lumb.kept.DE.gene.heatmap.pdf',width = 20,height=20)
Heatmap(tmp[flag,], col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = TRUE,show_heatmap_legend = FALSE,cluster_columns =  FALSE)
dev.off()

tmp           <- apply((median.tpm.matrix[her2.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
flag          <- which(rowSums(tmp) > 0)
pdf(file= '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/her2.kept.DE.gene.heatmap.pdf',width = 20,height=20)
Heatmap(tmp[flag,], col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = TRUE,show_heatmap_legend = FALSE,cluster_columns =  FALSE)
dev.off()

tmp           <- apply((median.tpm.matrix[basal.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(median.tpm.matrix)
flag          <- which(rowSums(tmp) > 0)
pdf(file = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/basal.kept.DE.gene.heatmap.pdf',width = 20,height=20)
Heatmap(tmp[flag,], col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = TRUE,show_heatmap_legend = FALSE,cluster_columns =  FALSE)
dev.off()




##############################################################################
#Figure S2e: propotion of genes belonging to different gene categories
##############################################################################

immune.gene.list <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/immune.gene.list.csv", stringsAsFactors=FALSE)$x %>% as.character()
lumb.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv",  stringsAsFactors=FALSE)$x
basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
her2.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv",  stringsAsFactors=FALSE)$x
lumb.dn.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.dn.csv",  stringsAsFactors=FALSE)$x
her2.dn.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.dn.csv",  stringsAsFactors=FALSE)$x
basal.dn.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.dn.csv", stringsAsFactors=FALSE)$x

lumb.gene  <- c(lumb.up.gene,lumb.dn.gene)
her2.gene  <- c(her2.up.gene,her2.dn.gene)
basal.gene <- c(basal.up.gene,basal.dn.gene)


lumb.df   <- rbind(data.frame(subtype='LuminalB',     y=setdiff(lumb.gene,immune.gene.list) %>% length, class='non.immune.gene'),data.frame(subtype='LuminalB',     y=intersect(lumb.gene,immune.gene.list) %>% length, class='immune.gene'))
lumb.df$y <- lumb.df$y/sum(lumb.df$y)
her2.df   <- rbind(data.frame(subtype='Her2-enriched',y=setdiff(her2.gene,immune.gene.list) %>% length, class='non.immune.gene'),data.frame(subtype='Her2-enriched',y=intersect(her2.gene,immune.gene.list) %>% length, class='immune.gene'))
her2.df$y <- her2.df$y/sum(her2.df$y)
basal.df  <- rbind(data.frame(subtype='Basal-like',   y=setdiff(basal.gene,immune.gene.list) %>% length,class='non.immune.gene'),data.frame(subtype='Basal-like',   y=intersect(basal.gene,immune.gene.list) %>% length,class='immune.gene'))
basal.df$y <- basal.df$y/sum(basal.df$y)

draw.df  <- rbind(lumb.df,her2.df,basal.df)
draw.df$subtype <- factor(draw.df$subtype,levels = c('Basal-like','Her2-enriched','LuminalB'))

ggplot(draw.df) + geom_bar(aes(x=subtype,y=y,fill=class),stat='identity',position=position_dodge()) + ggplot.style + theme(legend.position = "none") + scale_fill_manual(values=c('black','grey')) + ylab('Count')

ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/gene.class.bar.plot.pdf',width = 20,height=20)













################################ Trash code ################################ 
# res <- de.res.liver.vs.breast
# de.res.liver.vs.breast.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.01]
# de.res.liver.vs.breast.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.01]
# 
# 
# 
# res <- de.res.metastasis.liver.vs.breast.basal
# de.res.metastasis.liver.vs.breast.basal.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.01]
# de.res.metastasis.liver.vs.breast.basal.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.01]
# 
# 
# 
# res <- de.res.metastasis.liver.vs.breast.her2
# de.res.metastasis.liver.vs.breast.her2.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.01]
# de.res.metastasis.liver.vs.breast.her2.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.01]
# 
# res <- de.res.metastasis.liver.vs.breast.luma
# de.res.metastasis.liver.vs.breast.luma.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.01]
# de.res.metastasis.liver.vs.breast.luma.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.01]
# 
# 
# 
# g <- intersect(rownames(de.res.liver.vs.breast),rownames(de.res.metastasis.liver.vs.breast.her2))
# plot(x=de.res.liver.vs.breast[g,2],y=de.res.metastasis.liver.vs.breast.her2[g,2],xlim=c(-15,15),ylim=c(-15,15),xlab='liver.vs.breast',ylab='metastatic.vs.primary')
# abline(v=c(-1.5,1.5))
# abline(h=c(-1.5,1.5))
# lines(c(-15,15),c(-15,15))
# dd <- data.frame(x=de.res.liver.vs.breast[g,2],y=de.res.metastasis.liver.vs.breast.her2[g,2])
# rownames(dd) <- g
# View(dd[dd$x > 1.5 & dd$y < -1.5,])
# View(dd[dd$x < -1.5 & dd$y > 1.5,])
# 
# 
# 
# g <- intersect(rownames(de.res.liver.vs.breast),rownames(de.res.metastasis.liver.vs.breast.basal))
# plot(x=de.res.liver.vs.breast[g,2],y=de.res.metastasis.liver.vs.breast.basal[g,2],xlim=c(-15,15),ylim=c(-10,10))
# abline(v=c(-1.5,1.5))
# abline(h=c(-1.5,1.5))
# dd <- data.frame(x=de.res.liver.vs.breast[g,2],y=de.res.metastasis.liver.vs.breast.basal[g,2])
# rownames(dd) <- g
# View(dd[dd$x > 1.5 & dd$y < -1.5,])
# View(dd[dd$x < -1.5 & dd$y > 1.5,])
# 
# 
# 
# MET500.sample <- MET500.breast.cancer.polyA.LumB.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# TCGA.sample   <- TCGA.breast.cancer.polyA.LumB.sample
# 
# rs.MET500 <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# rs.TCGA <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# 
# #View(rs.MET500$correlation.matrix)
# s1 <-  rs.MET500$correlation.matrix[,'EFM192A_BREAST']
# 
# s2 <-  rs.TCGA$correlation.matrix[,'EFM192A_BREAST']
# 
# # plot(x=s1,y=MET500.log2.fpkm.matrix['ENSG00000153563',MET500.sample])
# # 
# # GATA2 <- 'ENSG00000179348'
# # GATAD2A <- 'ENSG00000167491'
# y <- c(MET500.log2.fpkm.matrix['ENSG00000136826',MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix['ENSG00000136826',TCGA.sample])
# x <- c(s1,s2)
# plot(x,y)
# 
# color.vec <- ifelse(names(y) %in% names(s1),'red','black')
# plot(x,y,col=color.vec,pch=19)
# 
# ADAT3 <-  'ENSG00000213638'
# 
# 
# ##############
# shuffle.gene.set <- foreach(x=xCell.corrected.gene.set) %do% {
#   sample(rownames(cancer.data.for.xCell),length(x)) 
#   
# }
# 
# 
# shuffle.cancer.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],
#                                            shuffle.gene.set, method = "ssgsea",
#                                    ssgsea.norm = FALSE)
# x             <- c(purity.TCGA,purity.MET500)
# 
# shuffle.cancer.de.df <- foreach(i =1:nrow(cancer.ssgsea.scores),.combine='rbind') %do% {
#   df                <- data.frame(x=x,y=shuffle.cancer.ssgsea.scores[i,])
#   loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#   adjusted.expr     <- loess.fit$residuals
#   p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
#   delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
#   data.frame(cancer.delta=delta,cancer.p.value=p.value)
#   
# }
# rownames(shuffle.cancer.de.df)  <- rownames(cancer.ssgsea.scores)
# #cancer.de.df$cancer.fdr <- p.adjust(cancer.de.df$cancer.p.value,method='fdr')
# 
# 
# 
# xCell.data$signatures[[1]]@geneIds
# 
# 
# 
# 
# data                  <- BRACA.log2.fpkm.matrix
# data                  <- data[rownames(data) %in% mapping.df$ENSEMBL,]
# rownames(data)        <- mapping.df[rownames(data),'SYMBOL']
# hehe <- GSVA::gsva(expr=data,
#                    xCell.corrected.gene.set, method = "ssgsea",
#                    ssgsea.norm = FALSE)
# 
# 
# 
# 
# 
# MET500.sample <- MET500.breast.cancer.polyA.Basal.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# TCGA.sample   <- TCGA.breast.cancer.polyA.Basal.sample
# 
# cancer.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],
#                                    xCell.data$signatures, method = "ssgsea",
#                                    ssgsea.norm = FALSE) # hmm, maybe here should set some cutoff, ssgsea score smaller than this cut-off truncate to zero. Let us check CCLE data
# 
# #xCell.score.matrix <- xCellAnalysis(cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],file.name = 'client-side/output/TME.R.output/basal.raw.score.txt' )
# 
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500 <- rs.MET500$correlation.matrix[,'HCC70_BREAST'] # use HCC70 cell line
# purity.TCGA   <- rs.TCGA$correlation.matrix[,'HCC70_BREAST']
# x             <- c(purity.TCGA,purity.MET500)
# 
# 
# cancer.geneset.da.df <- foreach(i =1:nrow(cancer.ssgsea.scores),.combine='rbind') %do% {
#   df                <- data.frame(x=x,y=cancer.ssgsea.scores[i,])
#   loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#   adjusted.expr     <- loess.fit$residuals
#   p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
#   delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
#   data.frame(cancer.delta=delta,cancer.p.value=p.value)
#   
# }
# rownames(cancer.geneset.da.df)  <- rownames(cancer.ssgsea.scores)
# cancer.geneset.da.df$cancer.fdr <- p.adjust(cancer.geneset.da.df$cancer.p.value,method='fdr')
# da.gene.set                     <- rownames(cancer.geneset.da.df)[cancer.geneset.da.df$cancer.fdr < 0.05]
# 
# Basal.cancer.geneset.da.df      <- cancer.geneset.da.df
# Basal.da.gene.set               <- da.gene.set

# ############################################################################################################
# # FigS1e: using liner regression to estimate abundance of hepatocytes
# #############################################################################################################
# 
# load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')
# sample.id <- 'SRR4306649'
# df <- hepatocyte.abundance.df[[sample.id]]
# intercept <- median(df$met.expr - df$liver.expr)
# ggplot(df,aes(x=liver.expr,y=met.expr)) + geom_point(size=2.5) + ggplot.style  + geom_abline(intercept = intercept,slope=1,colour='red',size=2) + xlim(c(-3,14)) + ylim(c(-3,14))+ geom_abline(intercept = 0,slope=1,linetype='dashed',size=2) + xlab('GTEx.liver') + ylab('MET500.sample')
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/estimate.hepatocyte.abundance.pdf',width = 20,height=20)


# my_palette <- colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev
# color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
#   scale = (length(lut)-1)/(max-min)
#   
#   dev.new(width=1.75, height=5)
#   plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
#   axis(2, ticks, las=1)
#   for (i in 1:(length(lut)-1)) {
#     y = (i-1)/scale + min
#     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
#   }
# }
# dev.off()
# pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/color.bar.pdf',width=20,height=20)
# color.bar(my_palette,min = -5,max=5,nticks = 30) # draw the color bar
# dev.off()




