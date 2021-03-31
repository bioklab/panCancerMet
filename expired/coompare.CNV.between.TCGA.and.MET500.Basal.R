require(ggplot2)
require(dplyr)
require(foreach)
require(dplyr)

load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/TCGA.breast.cancer.meta.R.output/TCGA.breast.cancer.meta.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')
load('client-side/output/organize.TCGA.and.MET500.breast.cancer.cnv.data.R.output/organize.TCGA.and.MET500.breast.cancer.cnv.data.RData')

MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
MET500.sample        <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Basal.sample

MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))

NOTCH.gene <- c('NRARP','MAML1','NOTCH1','MYC','CHAC1','MIB2','NOTCH3','BCL6','PTP4A3')

gene.vec <- intersect(rownames(MET500.breast.cancer.cnv.matrix), rownames(TCGA.breast.cancer.cnv.matrix))

MET500.vs.TCGA.cnv.comparision.df <- foreach(g = gene.vec,.combine='rbind') %do% {
    df1                  <- data.frame(cnv= MET500.breast.cancer.cnv.matrix[g,MET500.subject.id], site= rep('MET500',length(MET500.subject.id)))
    df2                  <- data.frame(cnv= TCGA.breast.cancer.cnv.matrix  [g,TCGA.subject.id],   site= rep('TCGA',length(TCGA.subject.id)))
    draw.df              <- rbind(df1,df2)
    draw.df$site         <- factor(draw.df$site,levels = c('TCGA','MET500'))
    p.value              <- wilcox.test(MET500.breast.cancer.cnv.matrix[g,MET500.subject.id],TCGA.breast.cancer.cnv.matrix[g,TCGA.subject.id])$p.value
    diff                 <- median(MET500.breast.cancer.cnv.matrix[g,MET500.subject.id]) - median(TCGA.breast.cancer.cnv.matrix[g,TCGA.subject.id])
    data.frame(p.value,diff,met=median(MET500.breast.cancer.cnv.matrix[g,MET500.subject.id]))
}
rownames(MET500.vs.TCGA.cnv.comparision.df) <- gene.vec
MET500.vs.TCGA.cnv.comparision.df           <- MET500.vs.TCGA.cnv.comparision.df[complete.cases(MET500.vs.TCGA.cnv.comparision.df),]
MET500.vs.TCGA.cnv.comparision.df$fdr       <- p.adjust(MET500.vs.TCGA.cnv.comparision.df$p.value,method='fdr')
NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df <- MET500.vs.TCGA.cnv.comparision.df[NOTCH.gene,]

save(file='client-side/output/coompare.CNV.between.TCGA.and.MET500.Basal.R.output/coompare.CNV.between.TCGA.and.MET500.Basal.RData',list=c('NOTCH.gene.MET500.vs.TCGA.cnv.comparision.df','MET500.vs.TCGA.cnv.comparision.df'))
