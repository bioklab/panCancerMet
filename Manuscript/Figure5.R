require(ggplot2)
source('client-side/code/Manuscript/ggplot.style.R')


##################################################################################################################################
#--------------------------- Figure 5 ---------------------------------------#
##################################################################################################################################

#################################################################
# (a): Table to show DE gene clusters
#################################################################


#################################################################
# (b): screen shot of UCSC genome browser to show the genomic amplification on chr19 (BRCA.Basal)
#################################################################


#################################################################
# (c): Survival analysis 19p13.12-amp vs normal
#################################################################
library(cgdsr)
library(survival)
library(survminer)

mycgds  <-  CGDS('http://www.cbioportal.org/',verbose=TRUE)
symbols <- c("CASP14", "ILVBL", "NOTCH3", "BRD4", "AKAP8L", "WIZ") 


mycancerstudy  <-  "brca_metabric" 
mycaselist     <- "brca_metabric_all"

myclinicaldata  <-  getClinicalData(mycgds,mycaselist)
myclinicaldata  <-  myclinicaldata[myclinicaldata$HER2_STATUS == "Negative" & myclinicaldata$ER_STATUS == "Negative" &  myclinicaldata$PR_STATUS == "Negative",]
myclinicaldata  <-  myclinicaldata[is.na(myclinicaldata$TUMOR_STAGE) == FALSE,]
myclinicaldata  <-  myclinicaldata[myclinicaldata$TUMOR_STAGE <= 2,]


mygeneticprofile <- 'brca_metabric_cna'
mrna             <- getProfileData(mycgds,c(as.character(symbols)) ,mygeneticprofile,mycaselist)
mrna             <- mrna[complete.cases(mrna),]
mrna.matrix      <- as.matrix(mrna)


expr_avg = apply(mrna.matrix, 1, function(x) median(x, na.rm = T))
expr_clinic <-  merge(x= data.frame(barcode = rownames(myclinicaldata), myclinicaldata), 
                      y= data.frame(barcode = rownames(mrna.matrix),    infect_expr = expr_avg  ), 
                      by = "barcode" 
)

expr_clinic  <- expr_clinic[,c('infect_expr','OS_MONTHS','OS_STATUS')]
expr_clinic  <- expr_clinic[complete.cases(expr_clinic),]


expr_clinic$group <- ifelse( expr_clinic$infect_expr >= 1, "GAIN", "OTHER")


expr_clinic$OS_STATUS_BIN = 1
expr_clinic$OS_STATUS_BIN[expr_clinic$OS_STATUS == "LIVING"] = 0

expr_clinic$group <- factor(expr_clinic$group,levels =  c('OTHER','GAIN'))

my.fit1 <- survfit(Surv(OS_MONTHS, OS_STATUS_BIN) ~ group, data =  expr_clinic) 
km.style <- theme_bw(base_size = 55) + theme(axis.title = element_text( size=30, face="bold"),
                                             axis.text  = element_text( size=55, face="bold"),
                                             plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                             axis.line.x = element_line(colour = "black",size = 3),
                                             axis.line.y = element_line(colour = "black",size = 3),
                                             plot.title = element_text(size =40, face = "bold")
) 
p <- ggsurvplot(my.fit1,data =expr_clinic,  pval = F,legend='none',ggtheme =km.style, xlab='',ylab='',palette = c("blue", "red"),size=4)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/19p13.12.KM.plot.pdf',width=20,height=15)







#################################################################
# (d):  KM plot of HES4 in basal-like breast cancer
#################################################################

load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')


TCGA.CDR.data <- read_excel("client-side/Data/TCGA-CDR/TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR") %>% as.data.frame
get.patient.id <- function(x) { 
    tmp <- strsplit(x,split = '-') %>% unlist
    paste(tmp[1:3],collapse = '-')
}

PRI.log2.tpm.matrix    <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]

g                   <- 'ENSG00000188290'
v                   <- PRI.log2.tpm.matrix[g,] 
high.group          <- sapply(names(v)[ v >= median(v)],get.patient.id)
low.group           <- sapply(names(v)[ v <  median(v)],get.patient.id)
high.group.df       <- TCGA.CDR.data[TCGA.CDR.data$bcr_patient_barcode %in% high.group,]
low.group.df        <- TCGA.CDR.data[TCGA.CDR.data$bcr_patient_barcode %in% low.group,]
high.group.df$group <- 'high'
low.group.df$group  <- 'low'
c.df                <- rbind(high.group.df,low.group.df)
c.df$group          <- factor(c.df$group,levels = c('low','high'))
c.df <- c.df[,c('PFI.time','PFI','group')]
c.df$PFI.time <- c.df$PFI.time / 30 # survival time in months
fit                 <- survfit(Surv(PFI.time,PFI) ~ group , data=c.df)
km.style <- theme_bw(base_size = 55) + theme(axis.title = element_text( size=30, face="bold"),
                                             axis.text  = element_text( size=70, face="bold"),
                                             plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                             axis.line.x = element_line(colour = "black",size = 3),
                                             axis.line.y = element_line(colour = "black",size = 3),
                                             plot.title = element_text(size =40, face = "bold")
) 
p <- ggsurvplot(fit=fit,data=c.df,legend='none',ggtheme =km.style, xlab='',ylab='',palette = c("blue", "red"),title= 'HES4',size=5) 
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/HES4.KM.plot.pdf',width=20,height=15)
sdf                 <- survdiff(Surv(PFI.time,PFI) ~  group, data = c.df)
p.val               <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)















##################################################################################################################################
#-------------- Figure SX ------------------------------------------------------------------------------
##################################################################################################################################

######################################################
# (a) (b): distribution of cluster size
######################################################
load('client-side/output/co.DE.cluster.R.output/co.DE.cluster.RData')

df      <- table(BRCA.Basal.fake.cluster.rs$gene.number) %>% as.data.frame
df$Var1 <- as.integer(df$Var1 %>% as.character())
df$freq <- df$Freq / sum(df$Freq)
#df      <- rbind(df,data.frame(Var1 = 5, Freq =0, freq = 0))
ggplot(df,aes(x=Var1,y=freq)) + geom_bar(stat="identity") + ggplot.style + xlab('cluster size') + ylab('Frequency') + ylim(0,1) + scale_x_continuous(breaks=2: max(df$Var1))
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/BRCA.Basal.cluster.size.bar.plot.pdf',width=20,height=15)
ggplot(df[2:nrow(df),],aes(x=Var1,y=freq)) + geom_bar(stat="identity") + ggplot.style + xlab('cluster size') + ylab('Frequency') +  scale_x_continuous(breaks=2: max(df$Var1))
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/BRCA.Basal.cluster.size.bar.plot.zoom.in.pdf',width=20,height=15)

df      <- table(BRCA.LumB.fake.cluster.rs$gene.number) %>% as.data.frame
df$Var1 <- as.integer(df$Var1 %>% as.character())
df$freq <- df$Freq / sum(df$Freq)
ggplot(df,aes(x=Var1,y=freq)) + geom_bar(stat="identity") + ggplot.style + xlab('cluster size') + ylab('Frequency') + ylim(0,1) + scale_x_continuous(breaks=2: max(df$Var1))
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/BRCA.LumB.cluster.size.bar.plot.pdf',width=20,height=15)
ggplot(df[2:nrow(df),],aes(x=Var1,y=freq)) + geom_bar(stat="identity") + ggplot.style + xlab('cluster size') + ylab('Frequency') +  scale_x_continuous(breaks=2: max(df$Var1))
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/BRCA.LumB.cluster.size.bar.plot.zoom.in.pdf',width=20,height=15)





######################################################
# (c) (d): volcano plot of DE analysis (HES4.high vs HES4.low)
######################################################
load('client-side/output/NOTCH3-HES4.R.output/NOTCH3-HES4.RData')
NOTCH3     <- 'ENSG00000074181'

NOTCH3.df <- TCGA.HES4.high.vs.low.res[rownames(TCGA.HES4.high.vs.low.res) == NOTCH3,]
ggplot() + geom_point(data=TCGA.HES4.high.vs.low.res,aes(x=log2FoldChange,y= -1 * log10(padj)),size=2.5) + ggplot.style + ylim(0,15) + geom_point(data=NOTCH3.df,aes(x=log2FoldChange,y= -1 * log10(padj)),size=7.5,color='red')
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/TCGA.HES4.high.vs.low.de.pdf',width = 20,height=20)

NOTCH3.df <- SRP157974.HES4.high.vs.low.res[rownames(SRP157974.HES4.high.vs.low.res) == NOTCH3,]
ggplot() + geom_point(data=SRP157974.HES4.high.vs.low.res,aes(x=log2FoldChange,y= -1 * log10(padj)),size=2.5) + ggplot.style + ylim(0,15) + geom_point(data=NOTCH3.df,aes(x=log2FoldChange,y= -1 * log10(padj)),size=7.5,color='red')
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/SRP157974.HES4.high.vs.low.de.pdf',width = 20,height=20)


# #################################################################
# #### Fig 4a Compare NOTCH3 CNV between primary and metestatic cancer
# #################################################################
# load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
# load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
# load('client-side/output/TCGA.breast.cancer.meta.R.output/TCGA.breast.cancer.meta.RData')
# load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')
# load('client-side/output/organize.TCGA.and.MET500.breast.cancer.cnv.data.R.output/organize.TCGA.and.MET500.breast.cancer.cnv.data.RData')
# 
# MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# MET500.sample        <- MET500.breast.cancer.polyA.Basal.sample
# MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Basal.sample
# 
# MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
# TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))
# 
# df1                  <- data.frame(cnv= MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id], site= rep('MET500',length(MET500.subject.id)))
# df2                  <- data.frame(cnv= TCGA.breast.cancer.cnv.matrix  ['NOTCH3',TCGA.subject.id],   site= rep('TCGA',length(TCGA.subject.id)))
# draw.df              <- rbind(df1,df2)
# draw.df$site         <- factor(draw.df$site,levels = c('TCGA','MET500'))
# 
# ggplot(draw.df) + geom_boxplot(aes(x=site,y=cnv),outlier.shape=NA,lwd=3) + ggplot.style + geom_jitter(aes(x=site,y=cnv),size=5.5) + xlab('')
# wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/NOTCH3.cnv.basal.pdf',width = 20,height=20)
# 
# 
# #################################################################
# #### Fig 4b  generated by UCSC genome browser
# #################################################################
# 
# 
# 
# #################################################################
# #### Fig 4c  DE analysis between HES4.high vs HES4.low
# #################################################################

# 
# 
# 
# 
# 
# ########################################
# #--------------Figure S6---------------#
# ########################################
# 
# #################################################################
# #### Fig S6a  cnv value of NOTCH genes in MET500
# #################################################################
# load('client-side/output/coompare.CNV.between.TCGA.and.MET500.Basal.R.output/coompare.CNV.between.TCGA.and.MET500.Basal.RData')
# NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df <- NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df[order(NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df$fdr),]
# sig.NOTCH.gene <- rownames(NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df)[NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df$fdr < 0.05]
# 
# draw.df <- foreach(g= sig.NOTCH.gene,.combine = 'rbind') %do% {
#   data.frame(cnv= MET500.breast.cancer.cnv.matrix[g,MET500.subject.id], gene=g)
#   
# }
# draw.df$gene <- factor(draw.df$gene,levels=sig.NOTCH.gene %>% as.character())
# ggplot(draw.df)  + geom_point(aes(x=gene,y=cnv),size=5.5) + coord_flip() + ggplot.style + geom_abline(intercept = 0,slope = 0,linetype=2,size=2)
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/significant.NOTCH.gene.cnv.in.MET500.pdf',width = 20,height=20)
# 
# #################################################################
# #### Fig S6b  CNV-expression correlation
# #################################################################
# NOTCH3    <- 'ENSG00000074181'
# df.TCGA   <- data.frame (expr= TCGA.breast.cancer.log2.tpm.matrix[NOTCH3,TCGA.subject.id],cnv=TCGA.breast.cancer.cnv.matrix  ['NOTCH3',TCGA.subject.id])
# df.MET500 <- data.frame (expr= MET500.log2.tpm.matrix[NOTCH3,MET500.sample],              cnv=MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id])
# ggplot(rbind(df.MET500,df.TCGA)) + geom_point(aes(x=cnv,y=expr),size=5.5) + ggplot.style
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/NOTCH3.cnv.expr.correlation.pdf',width = 20,height=20)
# 
# 
# 
# #################################################################
# #### Fig S6c  boxplot of NOTCH3 cnv in Her2 and LuminalB subtype
# #################################################################
# MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# MET500.sample        <- MET500.breast.cancer.polyA.LumB.sample
# MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample          <- pure.TCGA.breast.cancer.polyA.LumB.sample
# MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
# TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))
# df1                  <- data.frame(cnv= MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id], site= rep('MET500',length(MET500.subject.id)))
# df2                  <- data.frame(cnv= TCGA.breast.cancer.cnv.matrix  ['NOTCH3',TCGA.subject.id],   site= rep('TCGA',length(TCGA.subject.id)))
# draw.df              <- rbind(df1,df2)
# draw.df$site         <- factor(draw.df$site,levels = c('TCGA','MET500'))
# ggplot(draw.df) + geom_boxplot(aes(x=site,y=cnv),outlier.shape=NA,lwd=2) + ggplot.style + geom_jitter(aes(x=site,y=cnv),size=5.5) + xlab('')
# wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/NOTCH3.cnv.lumb.pdf',width = 20,height=20)
# 
# 
# MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# MET500.sample        <- MET500.breast.cancer.polyA.Her2.sample
# MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Her2.sample
# MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
# TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))
# df1                  <- data.frame(cnv= MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id], site= rep('MET500',length(MET500.subject.id)))
# df2                  <- data.frame(cnv= TCGA.breast.cancer.cnv.matrix  ['NOTCH3',TCGA.subject.id],   site= rep('TCGA',length(TCGA.subject.id)))
# draw.df              <- rbind(df1,df2)
# draw.df$site         <- factor(draw.df$site,levels = c('TCGA','MET500'))
# ggplot(draw.df) + geom_boxplot(aes(x=site,y=cnv),outlier.shape=NA,lwd=2) + ggplot.style + geom_jitter(aes(x=site,y=cnv),size=5.5) + xlab('')
# wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/NOTCH3.cnv.her2.pdf',width = 20,height=20)
# 
# 
# 
# 
# ########################################################################
# ### Fig S6d valcano plot of SRP157974
# ############################################################################



########################################################################
### Fig S6e chip-seq plot, generated by UCSC genome browser
############################################################################



########### Trash ###########

# #################################################################
# #### Fig S5b  boxplot of HES1 and HES4 expr in  Basal-like subtype
# #################################################################
# 
# load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
# load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
# load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
# load('server-side/RData//Breast Invasive Carcinoma.RData')
# load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')
# 
# 
# TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
# MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# MET500.sample  <- MET500.breast.cancer.polyA.Basal.sample
# MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample
# 
# HES1 <- 'ENSG00000114315'
# df1 <- data.frame(expr=MET500.log2.fpkm.matrix[HES1,MET500.sample],condition='MET500')
# df2 <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[HES1,TCGA.sample],condition='TCGA')
# df <- rbind(df1,df2)
# df$condition <- factor(df$condition,levels = c('TCGA','MET500'))
# ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=3) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') 
# wilcox.test(df1$expr,df2$expr)
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/HES1.DE.pdf',width = 20,height=20)
# 
# 
# HES4 <- 'ENSG00000188290'
# df1 <- data.frame(expr=MET500.log2.fpkm.matrix[HES4,MET500.sample],condition='MET500')
# df2 <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[HES4,TCGA.sample],condition='TCGA')
# df <- rbind(df1,df2)
# df$condition <- factor(df$condition,levels = c('TCGA','MET500'))
# ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=3) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') 
# wilcox.test(df1$expr,df2$expr)
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/HES4.DE.pdf',width = 20,height=20)
# 
# ########################################################################
# ### Fig S5d  scatter plot NOTCH3-HES4 across CCLE TNBC cell lines
# ############################################################################
# load('~/Project/Cancer2CellLine/client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
# load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')
# 
# NOTCH3     <- 'ENSG00000074181'
# HES4       <- 'ENSG00000188290'
# c.cell.line <- intersect(CCLE.breast.cancer.Basal.cell.line,colnames(CCLE.log2.rpkm.matrix))
# CCLE.log2.rpkm.matrix[c(NOTCH3,HES4),c.cell.line] %>% t %>% plot
# 
# draw.df <- data.frame(NOTCH3=CCLE.log2.rpkm.matrix[c(NOTCH3),c.cell.line],HES4= CCLE.log2.rpkm.matrix[c(HES4),c.cell.line])
# ggplot(draw.df,aes(x=NOTCH3,y=HES4)) + geom_point(size=5.5) + ggplot.style + ylab('HES4') + xlab('NOTCH3')
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/HES4.NOTC3.scatter.plot.pdf',width = 20,height=20)
# cor(draw.df$NOTCH3,draw.df$HES4,method='spearman')
# 
# 
# #################################################################
# #### Fig S5e  boxplot of HES4 brain metastsis
# #################################################################
# MET500.brain.sample <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'BRAIN'] 
# MET500.sample       <- c(MET500.breast.cancer.polyA.Basal.sample)
# MET500.sample       <- intersect(MET500.sample,MET500.brain.sample)
# 
# TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample
# df1            <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[HES4,TCGA.sample],condition='TCGA')
# df2            <- data.frame(expr=MET500.log2.fpkm.matrix[HES4,MET500.sample],          condition='MET500')
# draw.df <- rbind(df1,df2)
# ggplot(draw.df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=1.2) + ggplot.style + geom_jitter(aes(x=condition,y=expr),size=4.5) + xlab('') + ylim(0,9)
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section4/HES4.brain.metastasis.pdf',width = 20,height=20)
# wilcox.test(df1$expr,df2$expr)

