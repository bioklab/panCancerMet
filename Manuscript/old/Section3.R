require(ggplot2)
require(org.Hs.eg.db)
require(AnnotationDbi)
source('client-side/code/Manuscript/ggplot.style.R')
load('server-side/RData/Breast Invasive Carcinoma.RData')
load('client-side/output/TCGA.breast.cancer.meta.R.output/TCGA.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')


TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
PRRX1                                     <- 'ENSG00000116132'

################################################################################   
#----------------------------- Figure 3 ---------------------------------------#
################################################################################  


##############################################
## Fig 3a, PRRX1-common.dn.gene co-expr 
##############################################
load('client-side/output/PRRX1.R.output/PRRX1.RData')

draw.df <- rbind(Basal.cor.df,LumB.cor.df,Her2.cor.df)
draw.df$subtype <- factor(draw.df$subtype,levels = c('Basal-like','LuminalB','Her2-enriched'))
ggplot(draw.df,aes(x=subtype,y=cor.value,fill=gene.type)) + geom_violin(lwd=3)  + 
  scale_fill_manual(values = c('other.gene'='grey','ECM.gene'='blue')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('')  + ylim (-1.2,1.2)





#############################################
### Fig 3b, PRRX1 boxplot 
#############################################
MET500.sample  <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample
df1            <- data.frame(expr=MET500.log2.fpkm.matrix[PRRX1,MET500.sample],condition='MET500')
df2            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[PRRX1,TCGA.sample],condition='TCGA')
df             <- rbind(df1,df2)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
wilcox.test(df1$expr,df2$expr)
df.basal         <- df
df.basal$subtype <- 'Basal-like'


MET500.sample  <- MET500.breast.cancer.polyA.Her2.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Her2.sample
df1            <- data.frame(expr=MET500.log2.fpkm.matrix[PRRX1,MET500.sample],condition='MET500')
df2            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[PRRX1,TCGA.sample],condition='TCGA')
df             <- rbind(df1,df2)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
wilcox.test(df1$expr,df2$expr)
df.her2         <- df
df.her2$subtype <- 'Her2-enriched'



MET500.sample  <- MET500.breast.cancer.polyA.LumB.sample
MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample
df1            <- data.frame(expr=MET500.log2.fpkm.matrix[PRRX1,MET500.sample],condition='MET500')
df2            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[PRRX1,TCGA.sample],condition='TCGA')
df             <- rbind(df1,df2)
df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
wilcox.test(df1$expr,df2$expr)
df.lumb         <- df
df.lumb$subtype <- 'LuminalB'

draw.df         <- rbind(df.basal,df.her2,df.lumb)
draw.df$subtype <- factor(draw.df$subtype,levels = c('Basal-like','LuminalB','Her2-enriched'))
ggplot(draw.df,aes(x=subtype,y=expr,fill=condition)) +
  geom_boxplot(lwd=1.5,outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(),size=3.0) + 
  scale_fill_manual(values = c('MET500'='red','TCGA'='grey')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('') 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/PRRX1.boxplot.pdf',width = 20,height=20)










# basal.dn.gene.id  <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
# basal.dn.gene.id  <- setdiff(basal.dn.gene.id,PRRX1)
# cor.matrix        <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.Basal.sample] ,log2.fpkm.matrix[,pure.TCGA.breast.cancer.polyA.Basal.sample] %>% t,method='spearman')
# cor.vec           <- c(cor.matrix)
# names(cor.vec)    <- rownames(log2.fpkm.matrix)
# cor.vec           <- sort(cor.vec,decreasing = TRUE)
# cor.vec           <- cor.vec[is.na(cor.vec) == FALSE]
# cor.vec           <- cor.vec[names(cor.vec) != PRRX1]
# tmp               <- setdiff(names(cor.vec),basal.dn.gene.id)
# boxplot(cor.vec[tmp],    cor.vec[basal.dn.gene.id])
# wilcox.test(cor.vec[tmp],cor.vec[basal.dn.gene.id])
# draw.df          <- rbind(data.frame(value=cor.vec[tmp],type='other.gene'),data.frame(value=cor.vec[basal.dn.gene.id],type='dn.gene'))
# draw.df$subtype  <- 'Basal-like'
# basal.draw.df    <- draw.df
# #ggplot(draw.df) + geom_boxplot(aes(x=type,y=value),lwd=3,outlier.shape = NA)+ ggplot.style + geom_jitter(aes(x=type,y=value),size=2.5,width=0.05) + xlab('') + ylim(-0.8,1.2)
# #ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/basal.co.expr.boxplot.pdf',width = 20,height=20)
# x                     <- cor.vec[basal.dn.gene.id]
# basal.x               <- x[x > 0.6]
# 
# 
# 
# 
# her2.dn.gene.id  <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
# her2.dn.gene.id  <- setdiff(her2.dn.gene.id,PRRX1)
# cor.matrix       <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.Her2.sample] ,log2.fpkm.matrix[,pure.TCGA.breast.cancer.polyA.Her2.sample] %>% t,method='spearman')
# cor.vec          <- c(cor.matrix)
# names(cor.vec)   <- rownames(log2.fpkm.matrix)
# cor.vec          <- sort(cor.vec,decreasing = TRUE)
# cor.vec          <- cor.vec[is.na(cor.vec) == FALSE]
# cor.vec          <- cor.vec[names(cor.vec) != PRRX1]
# tmp              <- setdiff(names(cor.vec),her2.dn.gene.id)
# wilcox.test(cor.vec[tmp],cor.vec[her2.dn.gene.id])
# boxplot(cor.vec[tmp],    cor.vec[her2.dn.gene.id])
# draw.df          <- rbind(data.frame(value=cor.vec[tmp],type='other.gene'),data.frame(value=cor.vec[her2.dn.gene.id],type='dn.gene'))
# draw.df$subtype  <- 'Her2-enriched'
# her2.draw.df    <- draw.df
# #ggplot(draw.df) + geom_boxplot(aes(x=type,y=value),lwd=3,outlier.shape = NA)+ ggplot.style + geom_jitter(aes(x=type,y=value),size=2.5,width=0.05) + xlab('') + ylim(-0.8,1.2)
# #ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/her2.co.expr.boxplot.pdf',width = 20,height=20)
# x                     <- cor.vec[her2.dn.gene.id]
# her2.x                <- x[x > 0.6]
# 
# 
# 
# lumb.dn.gene.id  <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
# lumb.dn.gene.id  <- setdiff(lumb.dn.gene.id,PRRX1)
# cor.matrix       <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.LumB.sample] ,log2.fpkm.matrix[,pure.TCGA.breast.cancer.polyA.LumB.sample] %>% t,method='spearman')
# cor.vec          <- c(cor.matrix)
# names(cor.vec)   <- rownames(log2.fpkm.matrix)
# cor.vec          <- sort(cor.vec,decreasing = TRUE)
# cor.vec          <- cor.vec[is.na(cor.vec) == FALSE]
# cor.vec          <- cor.vec[names(cor.vec) != PRRX1]
# tmp              <- setdiff(names(cor.vec),lumb.dn.gene.id)
# boxplot(cor.vec[tmp],    cor.vec[lumb.dn.gene.id])
# wilcox.test(cor.vec[tmp],cor.vec[lumb.dn.gene.id])
# draw.df          <- rbind(data.frame(value=cor.vec[tmp],type='other.gene'),data.frame(value=cor.vec[lumb.dn.gene.id],type='dn.gene'))
# draw.df$subtype  <- 'LuminalB'
# lumb.draw.df    <- draw.df
# #ggplot(draw.df) + geom_boxplot(aes(x=type,y=value),lwd=3,outlier.shape = NA)+ ggplot.style + geom_jitter(aes(x=type,y=value),size=2.5,width=0.05) + xlab('') + ylim(-0.8,1.2)
# #ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/lumb.co.expr.boxplot.pdf',width = 20,height=20)
# x                     <- cor.vec[lumb.dn.gene.id]
# lumb.x                <- x[x > 0.6]
# 
# 
# draw.df <- rbind(basal.draw.df,her2.draw.df,lumb.draw.df)
# draw.df$subtype <- factor(draw.df$subtype,levels = c('Basal-like','LuminalB','Her2-enriched'))
# ggplot(draw.df,aes(x=subtype,y=value,fill=type)) + geom_violin(lwd=3)  + 
#   scale_fill_manual(values = c('other.gene'='grey','dn.gene'='blue')) + 
#   theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
#                                    axis.text.y  = element_text( size=55, face="bold"),
#                                    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
#                                    axis.line.x = element_line(colour = "black",size = 3),
#                                    axis.line.y = element_line(colour = "black",size = 3),
#                                    axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
#                                    legend.position= 'none')            + xlab('')  + ylim (-0.7,1.2)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/co.expr.violin.plot.pdf',width = 20,height=20)


################################################################################   
#----------------------------- Figure S4 ---------------------------------------#
################################################################################  


####################################################################    
### Fig S4a: boxplot of PRRX1 in dataset GSE58078 
####################################################################    
load('server-side/RData/BRACA_SRP043470.RData')
pri.expr         <- BRACA_SRP043470_log2.fpkm.matrix[PRRX1,c('SRR1427482','SRR1427483','SRR1427484')]
liver.met.expr   <- BRACA_SRP043470_log2.fpkm.matrix[PRRX1,c('SRR1427487','SRR1427488','SRR1427489')]

df1            <- data.frame(expr=pri.expr,condition='PRI')
df2            <- data.frame(expr=liver.met.expr,condition='MET')
df             <- rbind(df1,df2)
df$condition   <- factor(df$condition,levels = c('PRI','MET'))
ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=2)  + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') + ggplot.style
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/PRRX1.boxplot.SRP043470.pdf',width = 20,height=20)
wilcox.test(df1$expr,df2$expr,paired=TRUE)
t.test(df1$expr,df2$expr,paired=TRUE)

# p-value two high, maybe due to sample size
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4289456/


###############################
### Fig S4b: boxplot of PRRX1 in brain metastsis and TCGA  
###############################
MET500.brain.sample <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'BRAIN'] 
MET500.sample       <- c(MET500.breast.cancer.polyA.Basal.sample)
MET500.sample       <- intersect(MET500.sample,MET500.brain.sample)

TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample
df1            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[PRRX1,TCGA.sample],condition='TCGA')
df2            <- data.frame(expr=MET500.log2.fpkm.matrix[PRRX1,MET500.sample],          condition='MET500')
draw.df <- rbind(df1,df2)
ggplot(draw.df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=1.2) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=4.5) + xlab('') + ylim(0,9)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/PRRX1.brain.metastasis.pdf',width = 20,height=20)
wilcox.test(df1$expr,df2$expr)



##################################
### Fig S4c : BP dot plot for genes with high PRRX1 co-expression
##################################
compute.jaccard.index  <- function(x,y){
  ( intersect(x,y) %>% length )/( unique(c(x,y))   %>% length )
  
}
get.gene.list <- function(x) {
  strsplit(x,split='\\/')[[1]] %>% unlist    
}

enriched.term.filtering <- function(x) {
  result    <- x %>% dplyr::filter(Count >=5 & qvalue <= 0.2 & pvalue<=0.001) %>% dplyr::arrange(pvalue) 
  if(nrow(result) == 0){
    return(NA)    
  }
  gene.str  <- result$geneID %>% unique
  idx.vec   <- sapply(gene.str,function(x) which(result$geneID == x) %>% min  )
  result    <- result[idx.vec,]
  gene.list <- lapply(result$geneID,get.gene.list)
  sim.matrix <- matrix(data=0,nrow=nrow(result),ncol=nrow(result))
  for(i in 1:nrow(sim.matrix)){
    for(j in 1:nrow(sim.matrix)){
      sim.matrix[i,j] <- compute.jaccard.index(gene.list[[i]],gene.list[[j]])    
    }
  }
  diag(sim.matrix) <- 0
  
  
  bit.vec   <- rep(1,nrow(result))
  for(i in 1:length(bit.vec)) {
    if(bit.vec[i] == 1){
      cover.index <- which(sim.matrix[i,] >= 0.4)
      bit.vec[cover.index] <- 0
    }        
  }
  result[bit.vec == 1,]
  
}

tmp                   <- c(names(lumb.x),names(her2.x),names(basal.x)) %>% table %>% as.data.frame
c.gene                <- tmp$.[tmp$Freq >= 2] %>% as.character()
c.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = c.gene)
c.gene.BP             <- enrichGO(gene=c.gene.annotation.df$SYMBOL %>% as.character(),keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
c.gene.BP.filtered    <- enriched.term.filtering(c.gene.BP)
c.gene.BP.filtered$Direction <- 'dn'

GO.enrichment.dotplot <- function(BP.df) {
  BP.df             <- arrange(BP.df,-1 * log10(pvalue))  
  BP.df$Description <- paste(BP.df$Description,sprintf('(%d)',BP.df$Count), sep =' ' )
  BP.df$y           <- BP.df$Description %>% as.character()
  BP.df$Description <- factor(BP.df$Description,levels=BP.df$y)
  BP.df$size        <- 4
  p <- ggplot(BP.df) +
    geom_point(aes(x=-1 * log10(pvalue), y=Description,colour=Direction,size=size )) + scale_size(range = c(1,12)) +  # weired trick, but working!
    scale_colour_manual(values=c('up'='red', 'dn'='blue')) +
    theme_bw(base_size = 55) + theme(axis.title   = element_text( size=25, face="bold"),
                                     axis.text.y  = element_text( size=40, face="bold"),
                                     axis.text.x  = element_text(size=35, face="bold"),
                                     plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                     axis.line.x = element_line(colour = "black",size = 3),
                                     axis.line.y = element_line(colour = "black",size = 3),
                                     legend.position = 'none')  + xlab('-log10(pvalue)') + ylab('')  + xlim(c(0,18))
  p
}
p <- GO.enrichment.dotplot(c.gene.BP.filtered)
ggsave(plot = p,filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/c.gene.dot.plot.pdf',width = 30,height=20)




############################ Trash code ##################################3

# ####################################################################    
# ### Fig S4d: boxplot of TWIST1 in TCGA and MET500
# ####################################################################  
# TWIST1         <- 'ENSG00000122691'
# 
# MET500.sample  <- MET500.breast.cancer.polyA.Basal.sample
# MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample
# df1            <- data.frame(expr=MET500.log2.fpkm.matrix[TWIST1,MET500.sample],condition='MET500')
# df2            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[TWIST1,TCGA.sample],condition='TCGA')
# df             <- rbind(df1,df2)
# df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
# #ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=3)  + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') + ggplot.style
# #ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/basal.TWIST1.boxplot.pdf',width = 20,height=20)
# wilcox.test(df1$expr,df2$expr)
# df.basal <- df
# df.basal$subtype <- 'Basal-like'
# 
# 
# MET500.sample  <- MET500.breast.cancer.polyA.Her2.sample
# MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Her2.sample
# df1            <- data.frame(expr=MET500.log2.fpkm.matrix[TWIST1,MET500.sample],condition='MET500')
# df2            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[TWIST1,TCGA.sample],condition='TCGA')
# df             <- rbind(df1,df2)
# df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
# #ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=3)  + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') + ggplot.style
# #ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/her2.TWIST1.boxplot.pdf',width = 20,height=20)
# wilcox.test(df1$expr,df2$expr)
# df.her2 <- df
# df.her2$subtype <- 'Her2-enriched'
# 
# 
# 
# MET500.sample  <- MET500.breast.cancer.polyA.LumB.sample
# MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample
# df1            <- data.frame(expr=MET500.log2.fpkm.matrix[TWIST1,MET500.sample],condition='MET500')
# df2            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[TWIST1,TCGA.sample],condition='TCGA')
# df             <- rbind(df1,df2)
# df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
# #ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=3)  + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') + ggplot.style
# #ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/lumb.TWIST1.boxplot.pdf',width = 20,height=20)
# wilcox.test(df1$expr,df2$expr)
# df.lumb <- df
# df.lumb$subtype <- 'LuminalB'
# 
# ggplot(rbind(df.basal,df.her2,df.lumb),aes(x=subtype,y=expr,fill=condition)) +
#   geom_boxplot(lwd=1.5,outlier.shape = NA) + 
#   geom_point(position=position_jitterdodge(),size=3.0) + 
#   scale_fill_manual(values = c('MET500'='red','TCGA'='grey')) + 
#   theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
#                                    axis.text.y  = element_text( size=55, face="bold"),
#                                    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
#                                    axis.line.x = element_line(colour = "black",size = 3),
#                                    axis.line.y = element_line(colour = "black",size = 3),
#                                    axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
#                                    legend.position= 'none')            + xlab('') 
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/TWIST1.boxplot.pdf',width = 20,height=20)








