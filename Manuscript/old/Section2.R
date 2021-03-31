require(ggplot2)
require(dplyr)
library(clusterProfiler)
library(enrichplot)
source('client-side/code/Manuscript/ggplot.style.R')
load('client-side/output/analyze.DE.gene.R.output/analyze.DE.gene.RData')
load('client-side/output/analyze.DE.gene.R.output/common.DE.gene.RData')


##############################################################################
#---------------- Figure 2 --------------------------------------------------#
##############################################################################

#######################################################
# Fig2a: VennDiagram to show DE gene overlappings
######################################################
lumb.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.up.gene.immune.excluded.csv",  stringsAsFactors=FALSE)$x %>% unique
basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x %>% unique
her2.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.up.gene.immune.excluded.csv",  stringsAsFactors=FALSE)$x %>% unique
lumb.dn.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.dn.gene.immune.excluded.csv",  stringsAsFactors=FALSE)$x %>% unique
her2.dn.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.dn.gene.immune.excluded.csv",  stringsAsFactors=FALSE)$x %>% unique
basal.dn.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x %>% unique

require(gplots)
venn(list(lumb=lumb.up.gene,her2=her2.up.gene,basal=basal.up.gene)) # Figures used in the manuscript were made by PPT according to this
venn(list(lumb=lumb.dn.gene,her2=her2.dn.gene,basal=basal.dn.gene)) # Figures used in the manuscript were made by PPT according to this



#######################################################
# Fig2b: GO enrichment on common downregulated genes
######################################################

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
                                     axis.text.y  = element_text( size=35, face="bold"),
                                     axis.text.x  = element_text(size=35, face="bold"),
                                     plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                     axis.line.x = element_line(colour = "black",size = 3),
                                     axis.line.y = element_line(colour = "black",size = 3),
                                     legend.position = 'none')  + xlab('-log10(pvalue)') + ylab('')  
  p
}
c.dn.gene.BP.filtered$Direction <- 'dn'

c.dn.gene.BP.dotplot <- GO.enrichment.dotplot(c.dn.gene.BP.filtered)
ggsave(plot = c.dn.gene.BP.dotplot,filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/c.dn.gene.BP.dotplot.pdf',width = 40,height=20)





#######################################################
# Fig2c: DE-gene stage association analysis
######################################################
load('client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData')

basal.stage.association.rs$up.score <- scale(basal.stage.association.rs$up.score,scale = TRUE, center=TRUE)
basal.stage.association.rs$dn.score <- scale(basal.stage.association.rs$dn.score,scale = TRUE, center=TRUE)
her2.stage.association.rs$up.score  <- scale(her2.stage.association.rs$up.score,scale  = TRUE, center=TRUE)
her2.stage.association.rs$dn.score  <- scale(her2.stage.association.rs$dn.score,scale  = TRUE, center=TRUE)
lumb.stage.association.rs$up.score  <- scale(lumb.stage.association.rs$up.score,scale  = TRUE, center=TRUE)
lumb.stage.association.rs$dn.score  <- scale(lumb.stage.association.rs$dn.score,scale  = TRUE, center=TRUE)


basal.stage.association.rs$subtype <- 'Basal-like'
her2.stage.association.rs$subtype  <- 'Her2-enriched'
lumb.stage.association.rs$subtype  <- 'LuminalB'

pooled.rs         <- rbind(basal.stage.association.rs,her2.stage.association.rs,lumb.stage.association.rs)
pooled.rs$subtype <- factor(pooled.rs$subtype,levels = c('Basal-like','LuminalB','Her2-enriched'))
pooled.rs$stage <- factor(pooled.rs$stage,levels = c('stage I','stage II','stage III','stage IV'))


ggplot(pooled.rs,aes(x=subtype,y=up.score,fill=stage))  + geom_boxplot(lwd=1.5,outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(),size=3.0) + 
  scale_fill_manual(values = c('stage I'='#EF8A62','stage II'='#D1E5F0','stage III' = '#2166AC', 'stage IV' = 'grey')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('') 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/up.gene.ssgsea.pdf',width = 20,height=20)

ggplot(pooled.rs,aes(x=subtype,y=dn.score,fill=stage))  + geom_boxplot(lwd=1.5,outlier.shape = NA) + 
  geom_point(position=position_jitterdodge(),size=3.0) + 
  scale_fill_manual(values = c('stage I'='#EF8A62','stage II'='#D1E5F0','stage III' = '#2166AC', 'stage IV' = 'grey')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('') 

ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/dn.gene.ssgsea.pdf',width = 20,height=20)




wilcox.test(basal.stage.association.rs$up.score[basal.stage.association.rs$stage == 'stage I'],basal.stage.association.rs$up.score[basal.stage.association.rs$stage == 'stage II'])
wilcox.test(basal.stage.association.rs$up.score[basal.stage.association.rs$stage == 'stage I'],basal.stage.association.rs$up.score[basal.stage.association.rs$stage == 'stage III'])
wilcox.test(basal.stage.association.rs$up.score[basal.stage.association.rs$stage == 'stage II'],basal.stage.association.rs$up.score[basal.stage.association.rs$stage == 'stage III'])


wilcox.test(her2.stage.association.rs$up.score[her2.stage.association.rs$stage == 'stage I'], her2.stage.association.rs$up.score[her2.stage.association.rs$stage == 'stage II'])
wilcox.test(her2.stage.association.rs$up.score[her2.stage.association.rs$stage == 'stage I'], her2.stage.association.rs$up.score[her2.stage.association.rs$stage == 'stage III'])
wilcox.test(her2.stage.association.rs$up.score[her2.stage.association.rs$stage == 'stage II'],her2.stage.association.rs$up.score[her2.stage.association.rs$stage == 'stage III'])


wilcox.test(lumb.stage.association.rs$up.score[lumb.stage.association.rs$stage == 'stage I'],lumb.stage.association.rs$up.score[lumb.stage.association.rs$stage == 'stage II'])
wilcox.test(lumb.stage.association.rs$up.score[lumb.stage.association.rs$stage == 'stage I'],lumb.stage.association.rs$up.score[lumb.stage.association.rs$stage == 'stage III'])
wilcox.test(lumb.stage.association.rs$up.score[lumb.stage.association.rs$stage == 'stage II'],lumb.stage.association.rs$up.score[lumb.stage.association.rs$stage == 'stage III'])


wilcox.test(basal.stage.association.rs$dn.score[basal.stage.association.rs$stage == 'stage I'],basal.stage.association.rs$dn.score[basal.stage.association.rs$stage == 'stage II'])
wilcox.test(basal.stage.association.rs$dn.score[basal.stage.association.rs$stage == 'stage I'],basal.stage.association.rs$dn.score[basal.stage.association.rs$stage == 'stage III'])
wilcox.test(basal.stage.association.rs$dn.score[basal.stage.association.rs$stage == 'stage II'],basal.stage.association.rs$dn.score[basal.stage.association.rs$stage == 'stage III'])


wilcox.test(her2.stage.association.rs$dn.score[her2.stage.association.rs$stage == 'stage I'],her2.stage.association.rs$dn.score[her2.stage.association.rs$stage == 'stage II'])
wilcox.test(her2.stage.association.rs$dn.score[her2.stage.association.rs$stage == 'stage I'],her2.stage.association.rs$dn.score[her2.stage.association.rs$stage == 'stage III'])
wilcox.test(her2.stage.association.rs$dn.score[her2.stage.association.rs$stage == 'stage II'],her2.stage.association.rs$dn.score[her2.stage.association.rs$stage == 'stage III'])


wilcox.test(lumb.stage.association.rs$dn.score[lumb.stage.association.rs$stage == 'stage I'], lumb.stage.association.rs$dn.score[lumb.stage.association.rs$stage == 'stage II'])
wilcox.test(lumb.stage.association.rs$dn.score[lumb.stage.association.rs$stage == 'stage I'], lumb.stage.association.rs$dn.score[lumb.stage.association.rs$stage == 'stage III'])
wilcox.test(lumb.stage.association.rs$dn.score[lumb.stage.association.rs$stage == 'stage II'],lumb.stage.association.rs$dn.score[lumb.stage.association.rs$stage == 'stage III'])



#######################################################
# Fig2d: Volcano plot for identification of metatsis driver and suppresor genes, Basal-like as example
######################################################
load('client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData')

flag    <- basal.survival.rs.up$p.value < 0.05 & basal.survival.rs.up$effect.size > 0
ggplot(basal.survival.rs.up) + geom_point(aes(x=effect.size, y = -1 * log10(p.value)),size=5.5) + ggplot.style + geom_point(data= basal.survival.rs.up[flag,],aes(x=effect.size, y = -1 * log10(p.value)),size=5.5,colour='red') + geom_abline(intercept =  -1 * log10(0.05),slope = 0,size = 2,linetype=2) + ylim(0,3)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/basal.up.gene.driver.volcano.pdf',width = 20,height=20)

flag    <- basal.survival.rs.dn$p.value < 0.05 & basal.survival.rs.dn$effect.size < 0
ggplot(basal.survival.rs.dn) + geom_point(aes(x=effect.size, y = -1 * log10(p.value)),size=5.5) + ggplot.style + geom_point(data= basal.survival.rs.dn[flag,],aes(x=effect.size, y = -1 * log10(p.value)),size=5.5,colour='red') + geom_abline(intercept =  -1 * log10(0.05),slope = 0,size = 2,linetype=2) + ylim(0,3)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/basal.dn.gene.suppressor.volcano.pdf',width = 20,height=20)





######################################################
#------------------Figure S4 ------------------------#
######################################################




######################################################
#Figure S4a,b: volcanot plot for her2-enriched and luminalB subtype#
######################################################

flag    <- her2.survival.rs.up$p.value < 0.05 & her2.survival.rs.up$effect.size > 0
ggplot(her2.survival.rs.up) + geom_point(aes(x=effect.size, y = -1 * log10(p.value)),size=5.5) + ggplot.style + geom_point(data= her2.survival.rs.up[flag,],aes(x=effect.size, y = -1 * log10(p.value)),size=5.5,colour='red') + geom_abline(intercept =  -1 * log10(0.05),slope = 0,size = 2,linetype=2) + ylim(0,3)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/her2.up.gene.driver.volcano.pdf',width = 20,height=20)
flag     <- her2.survival.rs.dn$p.value < 0.05 & her2.survival.rs.dn$effect.size < 0
ggplot(her2.survival.rs.dn) + geom_point(aes(x=effect.size, y = -1 * log10(p.value)),size=5.5) + ggplot.style + geom_point(data= her2.survival.rs.dn[flag,],aes(x=effect.size, y = -1 * log10(p.value)),size=5.5,colour='red') + geom_abline(intercept =  -1 * log10(0.05),slope = 0,size = 2,linetype=2) + ylim(0,3)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/her2.dn.gene.suppressor.volcano.pdf',width = 20,height=20)



flag    <- lumb.survival.rs.up$p.value < 0.05 & lumb.survival.rs.up$effect.size > 0
ggplot(lumb.survival.rs.up) + geom_point(aes(x=effect.size, y = -1 * log10(p.value)),size=5.5) + ggplot.style + geom_point(data= lumb.survival.rs.up[flag,],aes(x=effect.size, y = -1 * log10(p.value)),size=5.5,colour='red') + geom_abline(intercept =  -1 * log10(0.05),slope = 0,size = 2,linetype=2) + ylim(0,3)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/lumb.up.gene.driver.volcano.pdf',width = 20,height=20)
flag     <- lumb.survival.rs.dn$p.value < 0.05 & lumb.survival.rs.dn$effect.size < 0
ggplot(lumb.survival.rs.dn) + geom_point(aes(x=effect.size, y = -1 * log10(p.value)),size=5.5) + ggplot.style + geom_point(data= lumb.survival.rs.dn[flag,],aes(x=effect.size, y = -1 * log10(p.value)),size=5.5,colour='red') + geom_abline(intercept =  -1 * log10(0.05),slope = 0,size = 2,linetype=2) + ylim(0,3)
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/lumb.dn.gene.suppressor.volcano.pdf',width = 20,height=20)



######################################################
#Figure S4c: KM plot for S100B and FAT #
######################################################
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/InSilicoCRISPR/client-side/output/organize.TCGA.clinical.data.R.output/organize.TCGA.clinical.data.RData')
load('client-side/output/tumor.purity.based.on.cell.line.R.output//tumor.purity.based.on.cell.line.RData')
require(survival)
require(survminer)


TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
clinical.data                             <- TCGA.clinical.data.list$BRCA

draw.KM.plot <- function(g) {
  expr           <- TCGA.breast.cancer.log2.fpkm.matrix[g,subtype.sample]
  patient.name   <- gsub(expr %>% names,pattern = '\\-01',replacement = '')
  df <- data.frame(expr=expr,
                   status=clinical.data[patient.name,'vital_status'], 
                   d1=clinical.data[patient.name,'days_to_death'] ,
                   d2=clinical.data[patient.name,'days_to_last_followup'],
                   age= clinical.data[patient.name,'years_to_birth'],
                   pathologic_stage= clinical.data[patient.name,'pathologic_stage']
  )
  df[is.na(df$d1),'d1'] <- 0
  df[is.na(df$d2),'d2'] <- 0
  df$time <- (df$d1 + df$d2) / 30
  df      <- df[complete.cases(df),]
  q                          <- quantile(df$expr)
  df$group                   <- 'median'
  df[df$expr > q[3],'group'] <- 'high'
  df[df$expr <= q[3],'group']<- 'low'
  df                         <- df[df$group!= 'median',]
  df$group                   <- factor(df$group,levels = c('low','high'))
  fit                        <- survfit(Surv(time, status) ~ group , data=df)
  km.style <- theme_bw(base_size = 55) + theme(axis.title = element_text( size=30, face="bold"),
                                               axis.text  = element_text( size=70, face="bold"),
                                               plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                               axis.line.x = element_line(colour = "black",size = 3),
                                               axis.line.y = element_line(colour = "black",size = 3),
                                               plot.title = element_text(size =40, face = "bold")
  ) 
  ggsurvplot(fit=fit,data=df,legend='none',ggtheme = km.style,xlab='',ylab='',palette = c("blue", "red"),xlim=c(0,350),title=g,size=5,break.x.by=100) 
}

basal.eg.dn.gene         <- 'ENSG00000160307'  # S100B
basal.eg.up.gene         <- 'ENSG00000086570'  # FAT2
subtype.sample           <- intersect(pure.TCGA.breast.cancer.polyA.Basal.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
basal.eg.up.gene.KM.plot <- draw.KM.plot(basal.eg.up.gene)
basal.eg.dn.gene.KM.plot <- draw.KM.plot(basal.eg.dn.gene)

pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/basal.eg.up.gene.KM.plot.pdf',width=20,height=15)
print(basal.eg.up.gene.KM.plot[[1]])
dev.off()
basal.survival.rs.up[basal.eg.up.gene,]


pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/basal.eg.dn.gene.KM.plot.pdf',width=20,height=15)
print(basal.eg.dn.gene.KM.plot[[1]])
dev.off()
basal.survival.rs.dn[basal.eg.dn.gene,]




#######################################################
# FigS3 a,b,c: Visulaization of GO enrichment results for Basal
######################################################
load('client-side/output/analyze.DE.gene.R.output/analyze.DE.gene.RData')

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
                                     axis.text.y  = element_text( size=35, face="bold"),
                                     axis.text.x  = element_text(size=35, face="bold"),
                                     plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                     axis.line.x = element_line(colour = "black",size = 3),
                                     axis.line.y = element_line(colour = "black",size = 3),
                                     legend.position = 'none')  + xlab('-log10(pvalue)') + ylab('')  
  p
}
basal.dn.BP.filtered$Direction <- 'dn'
basal.up.BP.filtered$Direction <- 'up'
basal.BP.dotplot <- GO.enrichment.dotplot(rbind(basal.dn.BP.filtered,basal.up.BP.filtered))
ggsave(plot = basal.BP.dotplot,filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/basal.BP.dotplot.pdf',width = 40,height=20)


lumb.dn.BP.filtered$Direction <- 'dn'
lumb.up.BP.filtered$Direction <- 'up'
lumb.BP.dotplot <- GO.enrichment.dotplot(rbind(lumb.dn.BP.filtered,lumb.up.BP.filtered))
ggsave(plot = lumb.BP.dotplot,filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/lumb.BP.dotplot.pdf',width = 40,height=20)

#her2.up.BP.filtered is NA
her2.up.BP.filtered$Direction <- 'up'
her2.dn.BP.filtered$Direction <- 'dn'
her2.BP.dotplot <- GO.enrichment.dotplot(rbind(her2.dn.BP.filtered,her2.up.BP.filtered))
ggsave(plot = her2.BP.dotplot,filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/her2.BP.dotplot.pdf',width = 40,height=20)







####################################################### Trash ####################################################### 
#######################################################
# FigS2a: Bar plot to show subtype-sepcificity of DE genes
######################################################
# up.gene <- c(lumb.up.gene,her2.up.gene,basal.up.gene)
# dn.gene <- c(lumb.dn.gene,her2.dn.gene,basal.dn.gene)
# 
# up.gene.freq.df <- table(up.gene) %>% as.data.frame
# dn.gene.freq.df <- table(dn.gene) %>% as.data.frame
# up.gene.freq.df$color <- 'up'
# dn.gene.freq.df$color <- 'dn'
# colnames(up.gene.freq.df)[1] <- 'gene'
# colnames(dn.gene.freq.df)[1] <- 'gene'
# 
# ggplot(rbind(up.gene.freq.df,dn.gene.freq.df)) + 
#   geom_bar( aes( x=factor(Freq),fill=color),position= 'dodge') + 
#   ggplot.style + 
#   xlab('Number of subtypes') + 
#   scale_fill_manual(values=c('up'='red','dn'='blue'))
# 
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/DE.gene.subtype.specificity.pdf',width = 20,height=20)

# subtype.sample          <- intersect(pure.TCGA.breast.cancer.polyA.Her2.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
# her2.eg.up.gene.KM.plot <- draw.KM.plot(her2.eg.up.gene)
# her2.eg.dn.gene.KM.plot <- draw.KM.plot(her2.eg.dn.gene)
# 
# pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/her2.eg.up.gene.KM.plot.pdf',width=20,height=15)
# print(her2.eg.up.gene.KM.plot[[1]])
# dev.off()
# her2.survival.rs.up[her2.eg.up.gene,]
# 
# pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/her2.eg.dn.gene.KM.plot.pdf',width=20,height=15)
# print(her2.eg.dn.gene.KM.plot[[1]])
# dev.off()
# her2.survival.rs.dn[her2.eg.dn.gene,]
# 
# 
# 
# subtype.sample          <- intersect(pure.TCGA.breast.cancer.polyA.LumB.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
# lumb.eg.up.gene.KM.plot <- draw.KM.plot(lumb.eg.up.gene)
# lumb.eg.dn.gene.KM.plot <- draw.KM.plot(lumb.eg.dn.gene)
# 
# pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/lumb.eg.up.gene.KM.plot.pdf',width=20,height=15)
# print(lumb.eg.up.gene.KM.plot[[1]])
# dev.off()
# lumb.survival.rs.up[lumb.eg.up.gene,]
# 
# 
# pdf(file ='~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/lumb.eg.dn.gene.KM.plot.pdf',width=20,height=15)
# print(lumb.eg.dn.gene.KM.plot[[1]])
# dev.off()
# lumb.survival.rs.dn[lumb.eg.dn.gene,]

# 
# 
# 
# 
# 
# 
# 
# tmp <- basal.survival.rs.up[basal.survival.rs.up$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size > 0) / nrow(tmp)
# basal.eg.up.gene <- 'ENSG00000103253'
# 
# 
# tmp <- her2.survival.rs.up[her2.survival.rs.up$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size > 0) / nrow(tmp)
# her2.eg.up.gene <- 'ENSG00000156453'
# 
# 
# tmp <- lumb.survival.rs.up[lumb.survival.rs.up$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size > 0) / nrow(tmp)
# lumb.eg.up.gene <- 'ENSG00000132613'
# 
# 
# tmp <- basal.survival.rs.dn[basal.survival.rs.dn$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size < 0) / nrow(tmp)
# basal.eg.dn.gene <- 'ENSG00000160307'
# 
# 
# tmp <- her2.survival.rs.dn[her2.survival.rs.dn$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size < 0) / nrow(tmp)
# her2.eg.dn.gene <- 'ENSG00000182326'
# 
# 
# tmp <- lumb.survival.rs.dn[lumb.survival.rs.dn$p.value < 0.05 ,]
# tmp <- tmp[order(tmp$p.value),]
# sum(tmp$effect.size < 0) / nrow(tmp)
# lumb.eg.dn.gene <- 'ENSG00000188001'

