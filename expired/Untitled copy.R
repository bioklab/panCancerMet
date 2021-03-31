
#Aim: DE analysis between primary and metastatic breast cancer (liver) in subtype-specific manner

require(DESeq2)
require(dplyr)
library (VennDiagram)
require(gplots)
require(foreach)
require(segmented)
source('client-side/code/util.R')

protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character


################################################################################################################
#Load reference tissue data
################################################################################################################
load('server-side//RData//Liver.RData')
female.sample                     <- sample.meta.df$sample.id[sample.meta.df$gender == 'Female'] %>% as.character()
Ref.liver.log2.read.count.matrix  <- log2.read.count.matrix[,female.sample]
Ref.liver.log2.fpkm.matrix        <- log2.fpkm.matrix[,female.sample]
liver.expressed.gene              <- get.expressed.gene(Ref.liver.log2.fpkm.matrix)


# load('server-side/RData/SRP068976_PairNor.RData')
# Ref.liver.log2.read.count.matrix  <- SRP068976_log2.read.count.matrix 
# Ref.liver.log2.fpkm.matrix        <- SRP068976_log2.fpkm.matrix 
# liver.expressed.gene              <- get.expressed.gene(Ref.liver.log2.fpkm.matrix)




#######################################################################################################################################################################################
###### Prepare the data 
#######################################################################################################################################################################################
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')

TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix
TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 

#######################################################################################################################################################################################
###### Function to perform DE analysis between tumor samples
#######################################################################################################################################################################################
perform.DE.analysis.between.metastatic.and.primary.cancer <- function(){
  g1             <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample])
  g2             <- get.expressed.gene(MET500.log2.fpkm.matrix[,MET500.sample])
  expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
  
  expr.matrix    <- cbind(MET500.log2.read.count.matrix[expressed.gene,MET500.sample],TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample])
  expr.matrix    <- 2^expr.matrix - 1
  
  purity.MET500  <- tumor.purity.based.on.cell.line.vec[MET500.sample]
  purity.TCGA    <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
  
  df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
  df$purity      <- c(purity.MET500,purity.TCGA)
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





gene.filtering.1 <- function(){
    
    padj.cut.off      <- 0.05
    immune.gene.list  <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/immune.gene.list.csv", stringsAsFactors=FALSE)$x
  
    # get the DE genes
    up.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange > 1  & deseq2.res$padj < padj.cut.off]
    dn.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange < -1 & deseq2.res$padj < padj.cut.off]
  
    # piecewise linear regression
    c.gene        <- intersect(rownames(deseq2.ref.res),rownames(deseq2.res))
    x             <- deseq2.ref.res[c.gene,'log2FoldChange']
    y             <- deseq2.res[c.gene,'log2FoldChange']
    lin.mod       <- lm(y~x)
    segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)    
    tmp           <- summary(segmented.mod)
    psi           <- tmp$psi[1,'Est.']
    ref.specific.gene  <- rownames(deseq2.ref.res)[deseq2.ref.res$log2FoldChange >= psi]
  
  #up.gene.due.to.ref <- intersect(up.gene,ref.specific.gene)
  #up.gene            <- setdiff(up.gene,ref.specific.gene)
  
    ######### gene filtering based on correlation with purity#########
    p              <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
    cor.vec        <- cor(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample] %>% t,p,method='spearman') %>% c
    names(cor.vec) <- rownames(TCGA.breast.cancer.log2.fpkm.matrix)
    cor.data       <- cor.vec[c(up.gene,dn.gene)]
    model          <- normalmixEM(cor.data)
    mu.vec         <- model$mu
    sigma.vec      <- model$sigma
    if(mu.vec[1] < mu.vec[2]){
       mu     <- mu.vec[1]
       sigma  <- sigma.vec[1]
    } else {
      mu      <- mu.vec[2]
      sigma   <- sigma.vec[2]
    }
    cor.cut.off  <- mu + 3.0 * sigma
    neg.cor.gene <- names(cor.data)[ cor.data <= cor.cut.off]
    ggplot.graph <- ggplot(data.frame(x=cor.data), aes(x=x)) + geom_histogram() + geom_density(aes(y=..density.. * 50)) + geom_vline(xintercept = cor.cut.off)
  

    ####### gene filtering based on expression level #############
    # p                       <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
    # q                       <- quantile(p,seq(0,1,0.1))
    # most.pure.TCGA.sample   <- names(p)[ p > q['80%'] ]
    TCGA.median             <- apply(TCGA.breast.cancer.log2.fpkm.matrix[, TCGA.sample], 1, median)
    MET500.median           <- apply(MET500.log2.fpkm.matrix[, MET500.sample], 1, median)
    median.max              <- apply(rbind(TCGA.median,MET500.median),2,max)
    x                       <- median.max[immune.gene.list]
    x                       <- x[is.na(x) == FALSE]
    q                       <- quantile(x)
    fpkm.cut.off            <- q['75%'] + 1.5 * IQR(q)
    low.expr.gene           <- names(median.max)[median.max < fpkm.cut.off]
    
    # p                       <- tumor.purity.based.on.cell.line.vec[MET500.sample]
  # q                       <- quantile(p,seq(0,1,0.1))
  # most.pure.MET500.sample <- names(p)[ p > q['80%'] ]
  # if(length(most.pure.MET500.sample) == 1){
  #     MET500.median           <- MET500.log2.fpkm.matrix[c(up.gene,dn.gene), most.pure.MET500.sample] %>% c
  # }else{
  # }
  
  

  ####### gene filtering based on wilcoxon rank test #######################
  wilcox.test.p.value.vec <- foreach(g=c(up.gene,dn.gene),.combine='c') %do% {
      wilcox.test(MET500.log2.fpkm.matrix[g,MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix[g,TCGA.sample])$p.value  
  }
  names(wilcox.test.p.value.vec) <- c(up.gene,dn.gene)
  not.robust.gene <- names(wilcox.test.p.value.vec) [ wilcox.test.p.value.vec > 0.05]
  
  
  uDE.gene <- c(ref.specific.gene,low.expr.gene,neg.cor.gene,not.robust.gene) %>% unique()
  up.gene  <- setdiff(up.gene,uDE.gene)
  dn.gene  <- setdiff(dn.gene,uDE.gene)
  
  # ######### gene filtering based on expression level #######################
  # g             <- c(up.gene,dn.gene)
  # MET500.m.expr <- apply(MET500.log2.fpkm.matrix[,MET500.sample],1,median)
  # TCGA.m.expr   <- apply(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],1,median)
  # m.expr.matrix <-  rbind(MET500.m.expr,TCGA.m.expr)
  # max.expr      <- apply(m.expr.matrix,2,max)
  # flag          <- MET500.m.expr[up.gene] > fpkm.cut.off
  # up.gene       <- up.gene[flag]
  # flag          <- TCGA.m.expr[dn.gene] > fpkm.cut.off
  # dn.gene       <- dn.gene[flag]
  
  list(up.gene=up.gene,dn.gene=dn.gene,g=ggplot.graph)
}





load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
#######################################################################################################################################################################################
###### DE analysis, Basal subtype
#######################################################################################################################################################################################
MET500.sample  <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample

deseq2.res            <- de.res.metastasis.liver.vs.breast.basal
deseq2.ref.res        <- de.res.liver.vs.breast.basal
Basal.de.gene.rs.1    <- gene.filtering.1()


up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = Basal.de.gene.rs.1$up.gene) 
GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 & GO.rs.1.up$pvalue < 0.01, ]


dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = Basal.de.gene.rs.1$dn.gene) 
GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 & GO.rs.1.dn$pvalue < 0.01, ]



cor.matrix <- cor(TCGA.breast.cancer.log2.fpkm.matrix[Basal.de.gene.rs.1$dn.gene,TCGA.sample] %>% t,method='spearman')
rs <- get.spectral.clustering.coordinates(A= cor.matrix)




#######################################################################################################################################################################################
###### DE analysis, Her2 subtype
#######################################################################################################################################################################################
MET500.sample  <- MET500.breast.cancer.polyA.Her2.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Her2.sample


deseq2.res            <- de.res.metastasis.liver.vs.breast.her2
deseq2.ref.res        <- de.res.liver.vs.breast.her2
Her2.de.gene.rs.1     <- gene.filtering.1()



up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = Her2.de.gene.rs.1$up.gene) 
GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 & GO.rs.1.up$pvalue < 0.01, ]


dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = Her2.de.gene.rs.1$dn.gene) 
GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 & GO.rs.1.dn$pvalue < 0.01, ]




#######################################################################################################################################################################################
###### DE analysis, LuminalB subtype
#######################################################################################################################################################################################
MET500.sample  <- c(MET500.breast.cancer.polyA.LumB.sample)
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample

deseq2.res            <- de.res.metastasis.liver.vs.breast.lumb
deseq2.ref.res        <- de.res.liver.vs.breast.lumb
LumB.de.gene.rs.1     <- gene.filtering.1()



up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = LumB.de.gene.rs.1$up.gene) 
GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 & GO.rs.1.up$pvalue < 0.01, ]


dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = LumB.de.gene.rs.1$dn.gene) 
GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 & GO.rs.1.dn$pvalue < 0.01, ]


#######################################################################################################################################################################################
###### DE analysis, LuminalA subtype. However, results are NOT included in the manuscript. Too few MET500 samples
#######################################################################################################################################################################################
v <- c(Basal.de.gene.rs.1$up.gene,Her2.de.gene.rs.1$up.gene,LumB.de.gene.rs.1$up.gene) %>% table %>% as.data.frame()
c.gene <- v$.[v$Freq >= 2] %>% as.character()
c.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = c.gene) 

