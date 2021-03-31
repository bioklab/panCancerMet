require(CePa)
require(dplyr)
require(ComplexHeatmap)
require(RColorBrewer)
require(circlize)
require(segmented)
require(survminer)
require(clusterProfiler)
source('client-side/code/Manuscript/ggplot.style.R')
load('client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData')

##############################################################################
#-------------------------- Figure 3 ----------------------------------------#
##############################################################################




##############################################################################
#(a) and (b):  Stage-expression association
##############################################################################

stage.df$cancer.type <- factor(stage.df$cancer.type,levels = c('BRCA.Basal','BRCA.LumB','BRCA.Her2','COAD'))
ggplot(stage.df,aes(x=cancer.type,y=up.score,fill= stage )) +
  geom_boxplot(lwd=1.5,outlier.shape = NA) + ylim(-5,5) + 
  geom_point(position=position_jitterdodge(),size=3.0) + 
  scale_fill_manual(values = c('Stage I'='#EF8A62','Stage II'='#D1E5F0','Stage III' = '#2166AC', 'Stage IV' = 'grey')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('') 

ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/stage.up.ssgsea.score.pdf',width = 20,height=20)


ggplot(stage.df,aes(x=cancer.type,y=dn.score,fill= stage )) +
  geom_boxplot(lwd=1.5,outlier.shape = NA) + ylim(-5,5) + 
  geom_point(position=position_jitterdodge(),size=3.0) + 
  scale_fill_manual(values = c('Stage I'='#EF8A62','Stage II'='#D1E5F0','Stage III' = '#2166AC', 'Stage IV' = 'grey')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('') 

ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/stage.dn.ssgsea.score.pdf',width = 20,height=20)








##############################################################################
#Figure 2c  Survival analysis
##############################################################################
load('client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData')
PRAD.survival.rs.up <- survival.df[survival.df$cancer.type == 'PRAD' & survival.df$gene.type =='up',]
flag                 <- PRAD.survival.rs.up$p.val < 0.05 & PRAD.survival.rs.up$hr > 0
ggplot(PRAD.survival.rs.up,aes(x=hr, y = -1 * log10(p.val))) + geom_point(size=5.5) + ggplot.style + geom_point(data= PRAD.survival.rs.up[flag,],aes(x=hr, y = -1 * log10(p.val)),size=8,colour='red') + geom_abline(intercept =  -1 * log10(0.05),slope = 0,size = 2,linetype=2) + ylim(0,11) + xlim(-3,3)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/PRAD.up.gene.survival.volcano.plot.pdf',width = 20,height=20)

PRAD.survival.rs.dn <- survival.df[survival.df$cancer.type == 'PRAD' & survival.df$gene.type =='dn',]
flag                 <- PRAD.survival.rs.dn$p.val < 0.05 & PRAD.survival.rs.dn$hr < 0
ggplot(PRAD.survival.rs.dn,aes(x=hr, y = -1 * log10(p.val))) + geom_point(size=5.5) + ggplot.style + geom_point(data= PRAD.survival.rs.dn[flag,],aes(x=hr, y = -1 * log10(p.val)),size=8,colour='red') + geom_abline(intercept =  -1 * log10(0.05),slope = 0,size = 2,linetype=2) + ylim(0,11)+ xlim(-3,3)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/PRAD.dn.gene.survival.volcano.plot.pdf',width = 20,height=20)



driver.gene.df    <- survival.df[survival.df$p.val < 0.05 & survival.df$hr > 0 & survival.df$gene.type == 'up',]
suppresor.gene.df <- survival.df[survival.df$p.val < 0.05 & survival.df$hr < 0 & survival.df$gene.type == 'dn',]
table(driver.gene.df$cancer.type)
table(suppresor.gene.df$cancer.type)

PRAD.driver.GO             <- enrichGO(gene=driver.gene.df$symbol[driver.gene.df$cancer.type == 'PRAD'] %>% unique ,         keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
PRAD.suppresor.GO          <- enrichGO(gene=suppresor.gene.df$symbol[suppresor.gene.df$cancer.type == 'PRAD'] %>% unique ,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
PRAD.driver.GO             <- PRAD.driver.GO[PRAD.driver.GO$Count >=5,]
PRAD.suppresor.GO          <- PRAD.suppresor.GO[PRAD.suppresor.GO$Count >=5,]


Basal.driver.GO             <- enrichGO(gene=driver.gene.df$symbol[driver.gene.df$cancer.type == 'BRCA.Basal'] %>% unique ,      keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
Basal.suppresor.GO          <- enrichGO(gene=suppresor.gene.df$symbol[suppresor.gene.df$cancer.type == 'BRCA.Basal'] %>% unique ,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
Basal.driver.GO             <- Basal.driver.GO[Basal.driver.GO$Count >=5,]
Basal.suppresor.GO          <- Basal.suppresor.GO[Basal.suppresor.GO$Count >=5,]




sig.gene.ratio.df <- foreach(cancer.type = c('BRCA.Basal','BRCA.LumB','BRCA.Her2','PRAD','COAD'),.combine='rbind') %do% {
    flag.up  <- survival.df$cancer.type == cancer.type & survival.df$gene.type =='up' 
    ratio.up <- sum(survival.df$p.val[flag.up] < 0.05) / sum(flag.up)
    flag.dn  <- survival.df$cancer.type == cancer.type & survival.df$gene.type =='dn' 
    ratio.dn <- sum(survival.df$p.val[flag.dn] < 0.05) / sum(flag.dn)
    data.frame(cancer.type = cancer.type,ratio.up = ratio.up,ratio.dn = ratio.dn)
}



ratio.df <- foreach(cancer.type = c('BRCA.Basal','BRCA.LumB','BRCA.Her2','PRAD','COAD'),.combine='rbind') %do% {
    flag.up <- survival.df$cancer.type == cancer.type & survival.df$gene.type =='up' 
    ratio.up <- sum(survival.df$hr[flag.up] > 0 & survival.df$p.val[flag.up] < 0.05) / sum(flag.up)
    flag.dn <- survival.df$cancer.type == cancer.type & survival.df$gene.type =='dn' 
    ratio.dn <- sum(survival.df$hr[flag.dn] < 0 & survival.df$p.val[flag.dn] < 0.05) / sum(flag.dn)
    data.frame(cancer.type = cancer.type,ratio.up = ratio.up,ratio.dn = ratio.dn)
    
}


TCGA.CDR.data <- read_excel("client-side/Data/TCGA-CDR/TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR") %>% as.data.frame
get.patient.id <- function(x) { 
    tmp <- strsplit(x,split = '-') %>% unlist
    paste(tmp[1:3],collapse = '-')
}

load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')
load('server-side/RData//Prostate Adenocarcinoma.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')
PRI.log2.tpm.matrix          <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]
draw.KM.plot <- function(g){
    v                   <- PRI.log2.tpm.matrix[g,] 
    high.group          <- sapply(names(v)[ v >= median(v)],get.patient.id)
    low.group           <- sapply(names(v)[ v <  median(v)],get.patient.id)
    high.group.df       <- TCGA.CDR.data[TCGA.CDR.data$bcr_patient_barcode %in% high.group,]
    low.group.df        <- TCGA.CDR.data[TCGA.CDR.data$bcr_patient_barcode %in% low.group,]
    high.group.df$group <- 'high'
    low.group.df$group  <- 'low'
    c.df                <- rbind(high.group.df,low.group.df)
    c.df$group          <- factor(c.df$group,levels = c('low','high'))
    fit                 <- survfit(Surv(PFI.time,PFI) ~  group, data = c.df)

    km.style <- theme_bw(base_size = 55) + theme(axis.title = element_text( size=30, face="bold"),
                                             axis.text  = element_text( size=70, face="bold"),
                                             plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                             axis.line.x = element_line(colour = "black",size = 3),
                                             axis.line.y = element_line(colour = "black",size = 3),
                                             plot.title = element_text(size =40, face = "bold")
) 
    ggsurvplot(fit=fit,data=c.df,legend='none',ggtheme = km.style,xlab='',ylab='',palette = c("blue", "red"),title=g,size=5,xlim=c(0,5500)) 
}


flag.up <- survival.df$cancer.type == 'PRAD' & survival.df$gene.type =='up' & survival.df$p.val < 0.05
View(survival.df[flag.up,])

g1 <- 'ENSG00000143153' # gene with most negative hr, ATP1B1
g2 <- 'ENSG00000076382' # gene with most positive hr, SPAG5

plt <- draw.KM.plot(g1)
pdf(file ='~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/PRAD.ATP1B1.KM.plot.pdf',width=20,height=15)
print(plt[[1]])
dev.flush()
dev.off()

plt <- draw.KM.plot(g2)
pdf(file ='~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/PRAD.SPAG5.KM.plot.pdf',width=20,height=15)
print(plt[[1]])
dev.flush()
dev.off()



flag.dn <- survival.df$cancer.type == 'PRAD' & survival.df$gene.type =='dn' & survival.df$p.val < 0.05
View(survival.df[flag.dn,])

g1 <- 'ENSG00000197894'  # gene with most negative hr, ADH5
g2 <- 'ENSG00000245970' # gene with most positive hr, AP003352.1, a lncRNA

plt <- draw.KM.plot(g1)
pdf(file ='~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/PRAD.ADH5.KM.plot.pdf',width=20,height=15)
print(plt[[1]])
dev.flush()
dev.off()

plt <- draw.KM.plot(g2)
pdf(file ='~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/PRAD.AP003352.1.KM.plot.pdf',width=20,height=15)
print(plt[[1]])
dev.flush()
dev.off()


####################################################### Trash ####################################################### 
# #######################################################
# # Fig 2a: motivation to use piecewise linear regression to remove liver speicfic genes (LuminalB subtype)
# ######################################################
# c.gene        <- intersect(rownames(de.res.liver.vs.breast.lumb),rownames(de.res.metastasis.liver.vs.breast.lumb))
# x             <- de.res.liver.vs.breast.lumb[c.gene,'log2FoldChange']
# y             <- de.res.metastasis.liver.vs.breast.lumb[c.gene,'log2FoldChange']
# lin.mod       <- lm(y~x)
# segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
# tmp           <- summary(segmented.mod)
# psi           <- tmp$psi[1,'Est.']
# df            <- data.frame(x=de.res.liver.vs.breast.lumb[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.lumb[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
# rownames(df)  <- c.gene
# ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
# ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/lumb.log2FC.MP.and.normal.pdf',width = 20,height=20)
# 
# 
# #######################################################
# # Fig 2b: motivation to use piecewise linear regression to remove liver speicfic genes (GSE58708 dataset)
# ######################################################
# load('client-side/output/validation.of.confounding.R.output/validation.of.confounding.RData')
# 
# c.gene        <- intersect(rownames(SRP043470.de.res.liver.vs.breast),rownames(SRP043470.de.res.metastasis.liver.vs.breast))
# x             <- SRP043470.de.res.liver.vs.breast[c.gene,'log2FoldChange']
# y             <- SRP043470.de.res.metastasis.liver.vs.breast[c.gene,'log2FoldChange']
# lin.mod       <- lm(y~x)
# segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
# tmp           <- summary(segmented.mod)
# psi           <- tmp$psi[1,'Est.']
# df            <- data.frame(x=SRP043470.de.res.liver.vs.breast[c.gene,'log2FoldChange'],y=SRP043470.de.res.metastasis.liver.vs.breast[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
# rownames(df)  <- c.gene
# ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
# ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/SRP043470.log2FC.MP.and.normal.pdf',width = 20,height=20)
# 
# 
# #######################################################
# # Fig 2c: DEBoost pipeline
# ######################################################
# 
# 
# ##############################################################################
# #---------------- Figure S2 --------------------------------------------------#
# ##############################################################################
# 
# #######################################################
# # Fig S2a: motivation to use piecewise linear regression to remove liver speicfic genes (Her2 subtype)
# ######################################################
# c.gene        <- intersect(rownames(de.res.liver.vs.breast.her2),rownames(de.res.metastasis.liver.vs.breast.her2))
# x             <- de.res.liver.vs.breast.her2[c.gene,'log2FoldChange']
# y             <- de.res.metastasis.liver.vs.breast.her2[c.gene,'log2FoldChange']
# lin.mod       <- lm(y~x)
# segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
# tmp           <- summary(segmented.mod)
# psi           <- tmp$psi[1,'Est.']
# df            <- data.frame(x=de.res.liver.vs.breast.her2[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.her2[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
# rownames(df)  <- c.gene
# ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/her2.log2FC.MP.and.normal.pdf',width = 20,height=20)
# 
# 
# #######################################################
# # Fig S2b: motivation to use piecewise linear regression to remove liver speicfic genes (Basal-like subtype)
# ######################################################
# c.gene        <- intersect(rownames(de.res.liver.vs.breast.basal),rownames(de.res.metastasis.liver.vs.breast.basal))
# x             <- de.res.liver.vs.breast.basal[c.gene,'log2FoldChange']
# y             <- de.res.metastasis.liver.vs.breast.basal[c.gene,'log2FoldChange']
# lin.mod       <- lm(y~x)
# segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
# tmp           <- summary(segmented.mod)
# psi           <- tmp$psi[1,'Est.']
# df            <- data.frame(x=de.res.liver.vs.breast.basal[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.basal[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
# rownames(df)  <- c.gene
# ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5) + ggplot.style + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_vline(xintercept =psi,linetype=2,size=2) + geom_line(aes(x=x,y=fitted),col='red',lwd=3.5)
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section1/basal.log2FC.MP.and.normal.pdf',width = 20,height=20)
# 
# 











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
# tmp <- PRAD.survival.rs.dn[PRAD.survival.rs.dn$p.value < 0.05 ,]
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



# ##############################################################################
# #Figure 2a  GO enrichment analysis results
# ##############################################################################
# 
# load('client-side/output/DE.gene.overview.R.output/DE.gene.overview.RData')
# most.sig.GO.rs.list <- foreach(GO.rs = GO.rs.list) %do% {
#   list(up = GO.rs$up[1,], dn = GO.rs$dn[1,])
# }
# names(most.sig.GO.rs.list) <- names(GO.rs.list)
# 
# 
# df <- foreach(item = most.sig.GO.rs.list,.combine='rbind') %do% {
#   item$up
# }
# df$cancer.type <- names(most.sig.GO.rs.list)
# up.df          <- df
# up.df$class    <- 'up'
# 
# 
# df <- foreach(item = most.sig.GO.rs.list,.combine='rbind') %do% {
#   item$dn
# }
# df$cancer.type <- names(most.sig.GO.rs.list)
# dn.df          <- df
# dn.df$class    <- 'dn'
# 
# 
# 
# draw.df <- rbind(up.df,dn.df)
# 
# 
# draw.df$cancer.type <- factor(draw.df$cancer.type, levels = c('BRCA.Basal','BRCA.LumB','BRCA.Her2','PRAD','COAD','PNET','SINET'))
# 
# 
# ggplot(draw.df,aes(x=cancer.type,y= -1 * log10(pvalue),fill=class))  + 
#   geom_bar(stat="identity",position="dodge",width= 0.7,show.legend = FALSE) + 
#   coord_flip() + theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
#                                                   axis.text  = element_text( size=25, face="bold"),
#                                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
#                                                   axis.line.x = element_line(colour = "black",size = 3),
#                                                   axis.line.y = element_line(colour = "black",size = 3))  + 
#   scale_fill_manual(values = c('up' = 'red', 'dn' = 'blue'))
# 
# ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section2/GO.bar.plot.pdf',width = 20,height=20)



