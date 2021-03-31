require(ggplot2)
require(CePa)
require(dplyr)
require(ComplexHeatmap)
require(RColorBrewer)
require(circlize)
require(segmented)
require(foreach)
source('client-side/code/Manuscript/ggplot.style.R')
load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')


################################################################################################################################################################      
#--------------------------------------- Figure 2 ---------------------------------------
################################################################################################################################################################     





################################################
# panel (a): expression pattern of the highly-upregulated (log2fc > 5) genes (derived from metastasis.vs.primary, LuminalB subtype)
################################################
GTEx.median.tpm.matrix <- read.gct(file='client-side/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
get.gene.id <- function(x) {
    strsplit(x=x,split='\\.') %>% unlist %>% head(1)
}
rownames(GTEx.median.tpm.matrix) <- sapply(rownames(GTEx.median.tpm.matrix),get.gene.id)

up.gene   <- rownames(BRCA.LumB.DE.rs$deseq2.M.vs.P.res)[BRCA.LumB.DE.rs$deseq2.M.vs.P.res$log2FoldChange > 5  & BRCA.LumB.DE.rs$deseq2.M.vs.P.res$padj < 0.05]
up.gene   <- intersect(up.gene,GTEx.median.tpm.matrix %>% rownames)
dev.off()
pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/BRCA.LumB.up.gene.expression.pattern.across.GTEx.tissues.pdf',width = 20,height=20)
tmp           <- apply((GTEx.median.tpm.matrix[up.gene,]+1) %>% log2,1,scale) %>% t
colnames(tmp) <- colnames(GTEx.median.tpm.matrix)
Heatmap(tmp, col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = FALSE,show_heatmap_legend = FALSE)
dev.off()

draw.color.bar <- function(value.range, color.range){
  f <- colorRamp2(value.range,color.range)
  x <- sapply(value.range, f)
  plot(c(0, 40*length(x)), c(0, 100), type = "n")
  i <- 0:(length(x)-1)
  rect(i*40, 0, (i+1)*40, 100, col=x)
}
pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/color.bar.pdf',width = 20,height=20)
draw.color.bar(value.range = seq(from=-5,to=5,by=0.5),
               color.range = colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev
)
dev.off()

tissue.enriched.gene.df <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/tissue.enriched.gene.csv", stringsAsFactors=FALSE)
liver.specific.gene <- tissue.enriched.gene.df$ensg_id[tissue.enriched.gene.df$enriched.tissue == 'liver'] %>% as.character()
(intersect(liver.specific.gene,up.gene) %>% length) / length(up.gene)

################################################################################################
# panel (b): an example of outlier SNRPG
################################################################################################
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')


PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.LumB.sample]
MET.log2.tpm.matrix       <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.LumB.sample]

# This is how I pick out the outlier gene 
p.value.vec               <- foreach(g= rownames(PRI.log2.tpm.matrix)) %do% {
    rs <- wilcox.test(PRI.log2.tpm.matrix[g,], MET.log2.tpm.matrix[g,])
    rs$p.value
}
names(p.value.vec) <- rownames(PRI.log2.tpm.matrix)

nde.gene <- names(p.value.vec)[p.value.vec > 0.05]

obj     <- BRCA.LumB.DE.rs$deseq2.M.vs.P.res
up.gene <- rownames(obj)[obj$log2FoldChange > 1 & obj$padj < 0.05]
x       <- intersect(up.gene,nde.gene)
obj[x,] %>% View



df             <- data.frame(condition=c(rep(x='MET500',times=length(pure.MET.breast.cancer.LumB.sample)), rep(x='TCGA',times=length(pure.PRI.breast.cancer.LumB.sample))))
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
f1             <- df$condition == 'MET500'
f2             <- df$condition == 'TCGA'
g              <- 'ENSG00000143977' # 
df$expr        <- c(MET.log2.tpm.matrix[g,],PRI.log2.tpm.matrix[g,])  
wilcox.test.p.value <- wilcox.test(df$expr[f1],df$expr[f2])$p.value
ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=2) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') 

ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/SNRPG.boxplot.pdf',width = 20,height=20)

BRCA.LumB.DE.rs$deseq2.M.vs.P.res[g,'pvalue']
BRCA.LumB.DE.rs$deseq2.M.vs.P.res[g,'padj']
wilcox.test.p.value

source('client-side/code/re.perform.DE.for.SNRPG.R') # re-perform DE anlaysis after removing outlier, then check the p-value
BRCA.LumB.DE.rs.without.outlier[g,'pvalue']
BRCA.LumB.DE.rs.without.outlier[g,'padj']




################################################################################################
# panel (c): DEBoost pipeline + log2FC scatter plot of BRCA.LumB
################################################################################################

load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
de.res.liver.vs.breast.lumb            <- BRCA.LumB.DE.rs$deseq2.R.vs.P.res
de.res.metastasis.liver.vs.breast.lumb <- BRCA.LumB.DE.rs$deseq2.M.vs.P.res
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
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/BRCA.LumB.log2FC.MP.and.normal.pdf',width = 20,height=20)




########################################    
# (d): proportion of DESeq2-called DE genes retained byDEBoost
#########################################   
DE.rs.list        <- list(BRCA.LumB.DE.rs,BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, PRAD.DE.rs,COAD.DE.rs,NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCA.LumB',      'BRCA.Basal',      'BRCA.Her2',     'PRAD',    'COAD',   'PNET',         'SINET')

df <- foreach(cancer.type = names(DE.rs.list),.combine='rbind') %do% {
  DE.rs         <- DE.rs.list[[cancer.type]]
  DE.rs.M.vs.P  <- DE.rs$deseq2.M.vs.P.res
  before.cnt    <- rownames(DE.rs.M.vs.P)[DE.rs.M.vs.P$log2FoldChange > 1 & DE.rs.M.vs.P$padj < 0.05] %>% length
  after.cnt     <- length(DE.rs$tumor.intrinsic.DE.gene.rs$up.gene) + length(DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
  #r.propotion   <- after.cnt / before.cnt
  #e.propotion   <- 1 - r.propotion
  r.propotion   <- after.cnt 
  e.propotion   <- before.cnt - after.cnt
  
  data.frame(cancer.type = cancer.type, gene.type = c('excluded','retained'),value=c(e.propotion,r.propotion))
}
ggplot(df,aes(x=cancer.type,y=value,fill = gene.type)) + geom_bar(position = 'stack',stat = "identity") + ggplot.style + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank()) + scale_fill_manual(values=c('retained' = '#F05C5E','excluded' = '#EFCE9E' ))
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/DEBoost.filtered.gene.proportion.pdf',width = 20,height=20)




############################################################################################      
#--------------------------------------- Figure S1 ---------------------------------------#
###########################################################################################    

#########
# (a): heatmap to show expression pattern (across GTEx tissues) of highly-upregulated genes 
##########
load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')

DE.rs.list        <- list(BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, PRAD.DE.rs,COAD.DE.rs,NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCA.Basal', 'BRCA.Her2','PRAD', 'COAD', 'PNET', 'SINET')

up.gene.df <- foreach(cancer.type=names(DE.rs.list),.combine='rbind') %do% {
    DE.rs     <- DE.rs.list[[cancer.type]]
    up.gene   <- rownames(DE.rs$deseq2.M.vs.P.res)[DE.rs$deseq2.M.vs.P.res$log2FoldChange > 5  & DE.rs$deseq2.M.vs.P.res$padj < 0.05]
    up.gene   <- intersect(up.gene,GTEx.median.tpm.matrix %>% rownames)
    data.frame(cancer.type = cancer.type,up.gene = up.gene)
}

tmp        <- apply((GTEx.median.tpm.matrix[up.gene.df$up.gene %>% as.character() ,]+1) %>% log2,1,scale) %>% t
cnt        <- apply(tmp,1, function(x) is.na(x) %>% sum)
idx        <- which(cnt > 0)
colnames(tmp) <- colnames(GTEx.median.tpm.matrix)
if(length(idx) > 0 ){
    tmp <- tmp[-1 * idx,]
}

row_anno <- rowAnnotation(cancer.type.block = anno_block())
pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/Six.cancer.type.pooled.up.gene.expression.pattern.across.GTEx.tissues.pdf',width = 20,height=20)
print(Heatmap(tmp,col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = FALSE,show_heatmap_legend = FALSE,cluster_rows = FALSE,right_annotation = row_anno,row_split = up.gene.df$cancer.type[-1 * idx] %>% as.integer(),row_gap = unit(10, "mm")))
dev.flush()
dev.off()


######################################     
#(b): boxplot to show CD8A differential expression 
######################################   
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')

PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.LumB.sample]
MET.log2.tpm.matrix       <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.LumB.sample]
df             <- data.frame(condition=c(rep(x='MET500',times=length(pure.MET.breast.cancer.LumB.sample)), rep(x='TCGA',times=length(pure.PRI.breast.cancer.LumB.sample))))
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
f1             <- df$condition == 'MET500'
f2             <- df$condition == 'TCGA'
g              <- 'ENSG00000153563' # CD8A
df$expr        <- c(MET.log2.tpm.matrix[g,],PRI.log2.tpm.matrix[g,])  
df$cancer.type <- 'BRCA.LumB'
df1 <- df
wilcox.test.p.value <- wilcox.test(df$expr[f1],df$expr[f2])$p.value
wilcox.test.p.value


# PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]
# MET.log2.tpm.matrix       <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.Basal.sample]
# df             <- data.frame(condition=c(rep(x='MET500',times=length(pure.MET.breast.cancer.Basal.sample)), rep(x='TCGA',times=length(pure.PRI.breast.cancer.Basal.sample))))
# df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
# f1             <- df$condition == 'MET500'
# f2             <- df$condition == 'TCGA'
# g              <- 'ENSG00000153563' # CD8A
# df$expr        <- c(MET.log2.tpm.matrix[g,],PRI.log2.tpm.matrix[g,])  
# df$cancer.type <- 'BRCA.Basal'
# df2 <- df
# wilcox.test.p.value <- wilcox.test(df$expr[f1],df$expr[f2])$p.value
# wilcox.test.p.value

ggplot(rbind(df1),aes(x=condition,y=expr)) + geom_boxplot(outlier.shape=NA,lwd=2,show.legend = FALSE) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('')  
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/CD8A.boxplot.pdf',width = 20,height=20)


################################################################################################################################################################      
#--------------------------------------- Figure S2 ---------------------------------------
################################################################################################################################################################     
for(cancer.type in names(DE.rs.list)) {
    DE.rs         <- DE.rs.list[[cancer.type]]
    DE.rs.R.vs.P  <- DE.rs$deseq2.R.vs.P.res
    DE.rs.M.vs.P  <- DE.rs$deseq2.M.vs.P.res
    c.gene        <- intersect(rownames(DE.rs.R.vs.P),rownames(DE.rs.M.vs.P))
    x             <- DE.rs.R.vs.P[c.gene,'log2FoldChange']
    y             <- DE.rs.M.vs.P[c.gene,'log2FoldChange']
    psi           <- -1
    while(psi < 0){ # well, for BRCA.Her2, this is necessary! psi = 6.877771
        lin.mod       <- lm(y~x)
        segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)    
        tmp           <- summary(segmented.mod)
        psi           <- tmp$psi[1,'Est.']
    }
    df            <- data.frame(x=DE.rs.R.vs.P[c.gene,'log2FoldChange'],y=DE.rs.M.vs.P[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
    rownames(df)  <- c.gene
    ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
    file.name <- sprintf('~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/%s.log2FC.MP.and.normal.pdf',cancer.type)
    ggsave(filename = file.name ,width = 20,height=20)
    print(sprintf('cancer.type = %s, psi = %f',cancer.type,psi))
}




################################################################################################################################################################      
#--------------------------------------- Figure S3 ---------------------------------------#
################################################################################################################################################################     








##################
# (b): heatmap to the expression pattern of up regulated genes retained by DEBoost
#################

load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')

DE.rs.list        <- list(BRCA.LumB.DE.rs,BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, PRAD.DE.rs,COAD.DE.rs,NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCAL.LumB',      'BRCA.Basal',     'BRCA.Her2',     'PRAD',    'COAD',    'PNET',        'SINET')

up.gene.df <- foreach(cancer.type=names(DE.rs.list),.combine='rbind') %do% {
  DE.rs     <- DE.rs.list[[cancer.type]]
  up.gene   <- DE.rs$tumor.intrinsic.DE.gene.rs$up.gene
  up.gene   <- intersect(up.gene,GTEx.median.tpm.matrix %>% rownames)
  data.frame(cancer.type = cancer.type,up.gene = up.gene)
}

tmp        <- apply((GTEx.median.tpm.matrix[up.gene.df$up.gene %>% as.character() ,]+1) %>% log2,1,scale) %>% t
cnt        <- apply(tmp,1, function(x) is.na(x) %>% sum)
idx        <- which(cnt > 0)
colnames(tmp) <- colnames(GTEx.median.tpm.matrix)
if(length(idx) > 0 ){
  tmp <- tmp[-1 * idx,]
}

row_anno <- rowAnnotation(cancer.type.block = anno_block())
pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/DEBoost.filtered.up.gene.expression.pattern.across.GTEx.tissues.pdf',width = 20,height=20)
print(Heatmap(tmp,col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = FALSE,show_heatmap_legend = FALSE,cluster_rows = FALSE,right_annotation = row_anno,row_split = up.gene.df$cancer.type[-1 * idx] %>% as.integer(),row_gap = unit(10, "mm")))
dev.flush()
dev.off()




########################################    
# (c): propotion of immune genes 
#########################################   
is.immune.gene <- function(x) {
  y    <- log2(x+1)    
  tmp  <- quantile(y)
  upper <- tmp[4] + 1.5* IQR(y)
  outlier.tissue <- names(y)[ y > upper]  
  immune.tissue    <- c('Spleen','Whole.Blood')
  if(sum(outlier.tissue %in% immune.tissue) ==2){
    return(TRUE)
  }else{
    return(FALSE)  
  }
}
flag               <- apply(GTEx.median.tpm.matrix,1,is.immune.gene)
immune.gene        <- rownames(GTEx.median.tpm.matrix)[flag]
names(immune.gene) <- NULL

DE.rs.list        <- list(BRCA.LumB.DE.rs,BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, PRAD.DE.rs,COAD.DE.rs,NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCA.LumB',      'BRCA.Basal',      'BRCA.Her2',     'PRAD',    'COAD',   'PNET',         'SINET')

df <- foreach(cancer.type = names(DE.rs.list),.combine='rbind') %do% {
    DE.rs   <- DE.rs.list[[cancer.type]]
    de.gene <- c(DE.rs$tumor.intrinsic.DE.gene.rs$up.gene, DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
    c.gene  <- intersect(immune.gene,de.gene)
    data.frame(cancer.type = cancer.type, proportion = length(c.gene) / length(de.gene))
}

ggplot(df) + geom_boxplot(aes(x='cancer',y=proportion),outlier.shape=NA,lwd=2) + ggplot.style + geom_jitter(aes(x='cancer',y=proportion),size=5.5) + xlab('')  + ylab('proportion ')
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/immune.gene.proportion.pdf',width = 20,height=20)

median(df$proportion)






























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
########################################    
#Fig S4a, boxplot of MYC expression
#########################################   
# PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]
# MET.log2.tpm.matrix       <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.Basal.sample]
# df             <- data.frame(condition=c(rep(x='MET500',times=length(pure.MET.breast.cancer.Basal.sample)), rep(x='TCGA',times=length(pure.PRI.breast.cancer.Basal.sample))))
# df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
# f1             <- df$condition == 'MET500'
# f2             <- df$condition == 'TCGA'
# g              <- 'ENSG00000136997' # MYC
# df$expr        <- c(MET.log2.tpm.matrix[g,],PRI.log2.tpm.matrix[g,])  
# wilcox.test.p.value <- wilcox.test(df$expr[f1],df$expr[f2])$p.value
# ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=2) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') 
# ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/BRCA.Basal.MYC.boxplot.pdf',width = 20,height=20)
# BRCA.Basal.DE.rs$deseq2.M.vs.P.res[g,'pvalue']
# wilcox.test.p.value
# 
# 
# 
