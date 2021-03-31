require(ggplot2)
require(dplyr)
source('client-side/code/Manuscript/ggplot.style.R')

load('client-side/output/pharmacological.interaction.R.output/pharmacological.interaction.RData')


####################################################################################
### Fig 5a, the pipeline
#####################################################################################


####################################################################################
### Fig 5b, scatterplot of RNAi and CRISPR residual score, BASAL cell line 
#####################################################################################

c.gene            <- intersect(rownames(BASAL.rs$CRISPR.result.df), rownames(BASAL.rs$RNAi.result.df))
draw.df           <- data.frame(x=BASAL.rs$CRISPR.result.df[c.gene,'sr'],y= BASAL.rs$RNAi.result.df[c.gene,'sr'])
rownames(draw.df) <- c.gene
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=8.5) + ggplot.style + ylab('RNAi') + xlab('CRISPR') + xlim(-10,10) + ylim(-10,10) 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/BASAL.screen.rs.pdf',width = 20,height=20)



####################################################################################
### Fig 5c, scatterplot of RNAi and CRISPR residual score, LYMPHOID 
#####################################################################################

c.gene            <- intersect(rownames(LYMPHOID.rs$CRISPR.result.df), rownames(LYMPHOID.rs$RNAi.result.df))
draw.df           <- data.frame(x=LYMPHOID.rs$CRISPR.result.df[c.gene,'sr'],y= LYMPHOID.rs$RNAi.result.df[c.gene,'sr'])
rownames(draw.df) <- c.gene
ggplot(draw.df,aes(x=x,y=y)) + geom_point(size=8.5) + ggplot.style + ylab('RNAi') + xlab('CRISPR') + xlim(-12,10) + ylim(-12,10) 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/LYMPHOID.screen.rs.pdf',width = 20,height=20)



####################################################################################
### Fig 5d, scatterplot of GDSC score 
#####################################################################################

ggplot(GDSC.IC50.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score)) + geom_point(size=8.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-9,9) + ylim(-9,9) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red') 
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/GDSC.screen.rs.pdf',width = 20,height=20)




##############################################################################################
### Fig S6a,b: scattler plot of median screen scroe in the two groups, BASAL breast cancer cell lines  
###############################################################################################

ggplot(BASAL.rs$CRISPR.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score))  + geom_point(size=3.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-6,6) + ylim(-6,6) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/BASAL.CRISPR.pdf',width = 20,height=20)

ggplot(BASAL.rs$RNAi.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score))  + geom_point(size=3.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-2,1.5) + ylim(-2,1.5) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/BASAL.RNAi.pdf',width = 20,height=20)



##############################################################################################
### Fig S6c,d: scattler plot of median screen scroe in the two groups, LYMPHOID cancer cell lines  
###############################################################################################

ggplot(LYMPHOID.rs$CRISPR.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score))  + geom_point(size=3.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-6,6) + ylim(-6,6) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/LYMPHOID.CRISPR.pdf',width = 20,height=20)


ggplot(LYMPHOID.rs$RNAi.result.df,aes(x=pan.cell.line.score,y=amp.cell.line.score))  + geom_point(size=3.5) + ggplot.style + ylab('amp.cell.line') + xlab('pan.cell.line') + xlim(-2,1.5) + ylim(-2,1.5) +  geom_abline(intercept = 0, slope = 1,linetype="dashed",size=2,colour='red')
ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/LYMPHOID.RNAi.pdf',width = 20,height=20)


############# Trash #######################
# load('client-side/output/pharmacological.interaction.R.output/pharmacological.interaction.RData')
# 
# ####################################################################################
# ## Fig 5a , show the negative correlation between notch3 pathway activity and TSSK6 depmap score 
# ####################################################################################
# cell.line <- intersect(names(NOTCH3.HES4.score),rownames(breast.cancer.crispr.screen.matrix))
# TSSK6.df <- data.frame(x=NOTCH3.HES4.score[cell.line],y=breast.cancer.crispr.screen.matrix[cell.line,'TSSK6'])
# ggplot(TSSK6.df,aes(x=x,y=y)) + geom_point(size=8.5) + ggplot.style + ylab('Depmap.score') + xlab('Pathway.activity') + xlim(0,70) 
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/TSSK6.crispr.screen.pdf',width = 20,height=20)
# cor(TSSK6.df$x,TSSK6.df$y,method='spearman')
# 
# 
# ######################################################
# ## Fig 5b , show the mean-variance relationship 
# ######################################################
# 
# ggplot(cp.with.replicates.df,aes(x=mean,y=sd)) + geom_point(size=8.5) + ggplot.style + ylab('Depmap.score') + xlab('Pathway.activity') + ylim(0,1)
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/sd.vs.mean.scatter.plot.pdf',width = 20,height=20)
# 
# 
# ################################################################################################
# ## Fig 5c， show the negative correlation between notch3 pathway activity and SN-38 depmap score #########
# ################################################################################################
# 
# df1                   <- data.frame(pert.id = primary.drug.screen.meta %>% rownames,   name = primary.drug.screen.meta$name,  dose = primary.drug.screen.meta$dose,  stringsAsFactors = FALSE)
# df2                   <- data.frame(pert.id = secondary.drug.screen.meta %>% rownames, name = secondary.drug.screen.meta$name,dose = secondary.drug.screen.meta$dose,stringsAsFactors = FALSE)
# screen.meta           <- rbind(df1,df2)
# screen.meta           <- screen.meta[complete.cases(screen.meta),]
# 
# SN38.pert.id          <- screen.meta$pert.id[screen.meta$name =='SN-38'] 
# SN38.cor.vec          <- drug.screen.notch3.hes4.cor.vec[SN38.pert.id] %>% sort 
# example.pert.id       <- names(SN38.cor.vec)[1]
# cell.line             <- intersect(names(NOTCH3.HES4.score),rownames(breast.cancer.secondary.drug.screen.matrix))
# SN38.df               <- data.frame(x=NOTCH3.HES4.score[cell.line],y=breast.cancer.secondary.drug.screen.matrix[cell.line,example.pert.id])
# 
# ggplot(SN38.df,aes(x=x,y=y)) + geom_point(size=8.5) + ggplot.style + ylab('Depmap.score') + xlab('Pathway.activity') + xlim(0,70) 
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/SN38.drug.screen.pdf',width = 20,height=20)
# cor(SN38.df$x,SN38.df$y,method='spearman')
# 
# ################################################################################################
# ## Fig 5d， show the negative correlation between notch3 pathway activity and exatecan.mesylate depmap score #########
# ################################################################################################
# 
# df1                   <- data.frame(pert.id = primary.drug.screen.meta %>% rownames,   name = primary.drug.screen.meta$name,  dose = primary.drug.screen.meta$dose,  stringsAsFactors = FALSE)
# df2                   <- data.frame(pert.id = secondary.drug.screen.meta %>% rownames, name = secondary.drug.screen.meta$name,dose = secondary.drug.screen.meta$dose,stringsAsFactors = FALSE)
# screen.meta           <- rbind(df1,df2)
# screen.meta           <- screen.meta[complete.cases(screen.meta),]
# 
# exatecan.mesylate.pert.id <- screen.meta$pert.id[screen.meta$name =='exatecan-mesylate'] 
# exatecan.mesylate.cor.vec <- drug.screen.notch3.hes4.cor.vec[exatecan.mesylate.pert.id] %>% sort 
# example.pert.id           <- names(exatecan.mesylate.cor.vec)[1]
# cell.line                 <- intersect(names(NOTCH3.HES4.score),rownames(breast.cancer.secondary.drug.screen.matrix))
# exatecan.mesylate.df      <- data.frame(x=NOTCH3.HES4.score[cell.line],y=breast.cancer.secondary.drug.screen.matrix[cell.line,example.pert.id])
# 
# ggplot(exatecan.mesylate.df,aes(x=x,y=y)) + geom_point(size=8.5) + ggplot.style + ylab('Depmap.score') + xlab('Pathway.activity') + xlim(0,70) 
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/exatecan.mesylate.drug.screen.pdf',width = 20,height=20)
# cor(exatecan.mesylate.df$x,exatecan.mesylate.df$y,method='spearman')
# 
# 
# #############################################
# ## Fig S6a depmap.crispr.score - NOTCH3.HES4.pathway.activity score correlation for SATL1 gene
# #############################################
# cell.line <- intersect(names(NOTCH3.HES4.score),rownames(breast.cancer.crispr.screen.matrix))
# SATL1.df <- data.frame(x=NOTCH3.HES4.score[cell.line],y=breast.cancer.crispr.screen.matrix[cell.line,'SATL1'])
# SATL1.df <- SATL1.df[complete.cases(SATL1.df),]
# ggplot(SATL1.df,aes(x=x,y=y)) + geom_point(size=8.5) + ggplot.style + ylab('Depmap.score') + xlab('Pathway.activity') + xlim(0,70) 
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/SATL1.crispr.screen.pdf',width = 20,height=20)
# cor(SATL1.df$x,SATL1.df$y,method='spearman')
# 
# 
# 
# 
# #############################################
# ## Fig S6b boxplot of TSSK6 expression in MET and TCGA
# ##############################################
# load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
# load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
# load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
# load('server-side/RData//Breast Invasive Carcinoma.RData')
# load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')
# 
# 
# TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
# MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# 
# 
# MET500.sample  <- MET500.breast.cancer.polyA.Basal.sample
# MET500.sample  <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample  <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample    <- pure.TCGA.breast.cancer.polyA.Basal.sample
# 
# TSSK6 <- 'ENSG00000178093'
# 
# df1 <- data.frame(expr=MET500.log2.fpkm.matrix[TSSK6,MET500.sample],condition='MET500')
# df2 <- data.frame(expr=TCGA.breast.cancer.log2.fpkm.matrix[TSSK6,TCGA.sample],condition='TCGA')
# df <- rbind(df1,df2)
# df$condition <- factor(df$condition,levels = c('TCGA','MET500'))
# ggplot(df) + geom_boxplot(aes(x=condition,y=expr),outlier.shape=NA,lwd=3)  + geom_jitter(aes(x=condition,y=expr),size=5.5) + xlab('') + ggplot.style + ylab('Expression')
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/TSSK6.expr.pdf',width = 20,height=20)
# 
# 
# 
# #############################################
# ## Fig S7 chenodeoxycholic.acid is an good example to show the in-consistency of the correlation values derived from different treatment 
# ##############################################
# 
# chenodeoxycholic.acid.pert.id    <- screen.meta$pert.id[screen.meta$name =='chenodeoxycholic-acid'] 
# chenodeoxycholic.acid.cor.vec    <- drug.screen.notch3.hes4.cor.vec[chenodeoxycholic.acid.pert.id] %>% sort
# id1 <- 'BRD-K18135438-001-14-2::2.5::MTS004'
# id2 <- 'BRD-K18135438-001-16-7::2.5::HTS'
# cell.line       <- intersect(names(NOTCH3.HES4.score),rownames(breast.cancer.secondary.drug.screen.matrix))
# chenodeoxycholic.acid.df         <- data.frame(x=NOTCH3.HES4.score[cell.line],
#                                         y1=breast.cancer.primary.drug.screen.matrix[cell.line,id1],
#                                         y2=breast.cancer.primary.drug.screen.matrix[cell.line,id2])
# 
# ggplot(chenodeoxycholic.acid.df,aes(x=x,y=y1)) + geom_point(size=8.5) + ggplot.style + ylab('Depmap.score') + xlab('Pathway.activity') 
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/chenodeoxycholic.acid.drug.screen.negative.pdf',width = 20,height=20)
# cor(chenodeoxycholic.acid.df$x,chenodeoxycholic.acid.df$y1,method='spearman',use='pairwise.complete.obs')
# 
# ggplot(chenodeoxycholic.acid.df,aes(x=x,y=y2)) + geom_point(size=8.5) + ggplot.style + ylab('Depmap.score') + xlab('Pathway.activity') 
# ggsave(filename = '~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/chenodeoxycholic.acid.drug.screen.positive.pdf',width = 20,height=20)
# cor(chenodeoxycholic.acid.df$x,chenodeoxycholic.acid.df$y2,method='spearman',use='pairwise.complete.obs')
# 
# 
# #############################################
# ######### Fig S8  SN-38 depmap.score - notch3.pathway correlation in multiple treatment
# #############################################
# 
# SN38.acid.pert.id    <- screen.meta$pert.id[screen.meta$name =='SN-38']
# cell.line            <- intersect(names(NOTCH3.HES4.score),rownames(breast.cancer.secondary.drug.screen.matrix))
# SN38.cor.vec         <- drug.screen.notch3.hes4.cor.vec[SN38.acid.pert.id] %>% sort
# SN38.cor.vec         <- SN38.cor.vec[2:length(SN38.cor.vec)]
# 
# require(foreach)
# fig.list <- foreach(pert.id = names(SN38.cor.vec)) %do% {
#   if(pert.id %in% colnames(breast.cancer.secondary.drug.screen.matrix)){
#       df         <- data.frame(x=NOTCH3.HES4.score[cell.line],
#                                y=breast.cancer.secondary.drug.screen.matrix[cell.line,pert.id]
#                               )
#   }else{
#       df         <- data.frame(x=NOTCH3.HES4.score[cell.line],
#                                y=breast.cancer.primary.drug.screen.matrix[cell.line,pert.id]
#                                )
#   }
#   ggplot(df,aes(x=x,y=y)) + geom_point(size=8.5) + ggplot.style + ylab('Depmap.score') + xlab('Pathway.activity') 
#   
# }
# for(i in 1:length(SN38.cor.vec)) {
#     ggsave(plot=fig.list[[i]],filename = sprintf('~/OneDrive/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section5/SN38.rep%d.drug.screen.positive.pdf',i),width = 20,height=20)
# }




