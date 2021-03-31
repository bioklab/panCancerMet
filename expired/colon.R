require(DESeq2)
require(dplyr)
require(foreach)
require(segmented)
source('client-side/code/util.R')
source('client-side/code/Manuscript/ggplot.style.R')

protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character


################################################################################################################
#Load reference tissue data
################################################################################################################
load('server-side//RData//Liver.RData')
Ref.liver.log2.fpkm.matrix        <- log2.fpkm.matrix
Ref.liver.log2.read.count.matrix  <- log2.read.count.matrix
liver.expressed.gene              <- get.expressed.gene(Ref.liver.log2.fpkm.matrix)



load('server-side//RData//COLORECTAL_SRP029880.RData')

TCGA.breast.cancer.log2.fpkm.matrix <- COLORECTAL_SRP029880_log2.fpkm.matrix
MET500.log2.fpkm.matrix             <- COLORECTAL_SRP029880_log2.fpkm.matrix

TCGA.breast.cancer.log2.read.count.matrix <- COLORECTAL_SRP029880_log2.read.count.matrix
MET500.log2.read.count.matrix             <- COLORECTAL_SRP029880_log2.read.count.matrix

TCGA.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'primary colorectal cancer'] %>% as.character()
MET500.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'metastatic colorectal cancer to the liver'] %>% as.character()


#######################################################################################################################################################################################
###### Function to perform DE analysis between tumor samples
#######################################################################################################################################################################################
perform.DE.analysis.between.metastatic.and.primary.cancer <- function(){
  g1             <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample])
  g2             <- get.expressed.gene(MET500.log2.fpkm.matrix[,MET500.sample])
  expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
  
  expr.matrix    <- cbind(MET500.log2.read.count.matrix[expressed.gene,MET500.sample],TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample])
  expr.matrix    <- 2^expr.matrix - 1
  

  df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
  df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
  
  dds           <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
                                          colData = df,
                                          design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds,contrast = c('condition','MET500','TCGA')) %>% as.data.frame
  res <- res[order(res$pvalue),]
  res <- res[complete.cases(res),]     
  res
}

perform.DE.analysis.between.ref.tissue.and.primary <- function(){
  g1             <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample])
  g2             <- get.expressed.gene(Ref.liver.log2.fpkm.matrix)
  expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
  
  expr.matrix    <- cbind(Ref.liver.log2.read.count.matrix[expressed.gene,],TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample])
  expr.matrix    <- 2^expr.matrix - 1
  
  df             <- data.frame(condition=c(rep(x='ref.tissue',times=ncol(Ref.liver.log2.read.count.matrix)), rep(x='TCGA',times=length(TCGA.sample))))
  df$condition   <- factor(df$condition,levels = c('TCGA','ref.tissue'))
  
  dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
                                           colData = df,
                                           design = ~ condition )
  dds <- DESeq(dds)
  res <- results(dds,contrast = c('condition','ref.tissue','TCGA')) %>% as.data.frame
  res <- res[order(res$pvalue),]
  res <- res[complete.cases(res),]     
  res
}



gene.filtering <- function(){
  c.gene        <- intersect(rownames(deseq2.ref.res),rownames(deseq2.res))
  x             <- deseq2.ref.res[c.gene,'log2FoldChange']
  y             <- deseq2.res[c.gene,'log2FoldChange']
  lin.mod       <- lm(y~x)
  segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)    
  tmp           <- summary(segmented.mod)
  psi           <- tmp$psi[1,'Est.']
  
  #ref.up.gene <- rownames(deseq2.ref.res)[deseq2.ref.res$log2FoldChange >  1  & deseq2.ref.res$padj < cut.off]
  #ref.dn.gene <- rownames(deseq2.ref.res)[deseq2.ref.res$log2FoldChange < -1  & deseq2.ref.res$padj < cut.off]
  
  up.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange > 1  & deseq2.res$padj < cut.off]
  dn.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange < -1 & deseq2.res$padj < cut.off]
  
  ref.specific.gene  <- rownames(deseq2.ref.res)[deseq2.ref.res$log2FoldChange >= psi]
  
  up.gene.due.to.ref <- intersect(up.gene,ref.specific.gene)
  up.gene            <- setdiff(up.gene,ref.specific.gene)
  
  ####### gene filtering based on wilcoxon rank test #######################
  wilcox.test.p.value.vec <- foreach(g=c(up.gene,dn.gene),.combine='c') %do% {
    wilcox.test(MET500.log2.fpkm.matrix[g,MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix[g,TCGA.sample])$p.value  
  }
  names(wilcox.test.p.value.vec) <- c(up.gene,dn.gene)
  up.gene <- up.gene[ wilcox.test.p.value.vec[up.gene] < 0.05 ]
  dn.gene <- dn.gene[ wilcox.test.p.value.vec[dn.gene] < 0.05 ]
  
  
  
  ######### gene filtering based on expression level #######################
  MET500.m.expr <- apply(MET500.log2.fpkm.matrix[,MET500.sample],1,median)
  TCGA.m.expr   <- apply(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],1,median)
  flag          <- MET500.m.expr[up.gene] > fpkm.cut.off
  up.gene       <- up.gene[flag]
  flag          <- TCGA.m.expr[dn.gene] > fpkm.cut.off
  dn.gene       <- dn.gene[flag]
  list(up.gene=up.gene,dn.gene=dn.gene)
}


fpkm.cut.off <- log2(10+1)
cut.off      <- 0.05

deseq2.res            <- perform.DE.analysis.between.metastatic.and.primary.cancer()
deseq2.ref.res        <- perform.DE.analysis.between.ref.tissue.and.primary()

require(ggplot2)
c.gene        <- intersect(rownames(deseq2.res),rownames(deseq2.ref.res))
x             <- deseq2.ref.res[c.gene,'log2FoldChange']
y             <- deseq2.res[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=deseq2.ref.res[c.gene,'log2FoldChange'],y=deseq2.res[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
df$residual   <- df$y - df$fitted
tissue.specific.gene.df <- df[df$x > psi,]

FOXA2 <- 'ENSG00000125798'
HNF1A <- 'ENSG00000135100'

tissue.specific.gene.df[FOXA2,]

tissue.specific.gene.df[HNF1A,]





de.gene.rs            <- gene.filtering()
basal.up.gene         <- de.gene.rs$up.gene
basal.dn.gene         <- de.gene.rs$dn.gene
de.res.metastasis.liver.vs.breast.basal         <- deseq2.res
de.res.liver.vs.breast.basal                    <- deseq2.ref.res
