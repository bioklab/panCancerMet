require(plyr)
require(dplyr)
require(foreach)
library(readxl)
library(writexl)

load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')

DE.rs.list        <- list(BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, BRCA.LumB.DE.rs, PRAD.DE.rs,COAD.DE.rs,NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCA.Basal', 'BRCA.Her2', 'BRCA.LumB','PRAD', 'COAD', 'PNET', 'SINET')




#Table S1
pooled.up.gene.df <- foreach(cancer.type = names(DE.rs.list),.combine='rbind') %do% {
    o          <- DE.rs.list[[cancer.type]] 
    up.gene    <- o$tumor.intrinsic.DE.gene.rs$up.gene      
    up.gene.df <- o$deseq2.M.vs.P.res[up.gene,]
    up.gene.df$cancer.type <- cancer.type
    up.gene.df$gene.id     <- up.gene
    up.gene.df
}
pooled.up.gene.df$direction <- 'upregulated'

pooled.dn.gene.df <- foreach(cancer.type = names(DE.rs.list),.combine='rbind') %do% {
  o          <- DE.rs.list[[cancer.type]] 
  dn.gene    <- o$tumor.intrinsic.DE.gene.rs$dn.gene      
  dn.gene.df <- o$deseq2.M.vs.P.res[dn.gene,]
  dn.gene.df$cancer.type <- cancer.type
  dn.gene.df$gene.id     <- dn.gene
  dn.gene.df
}
pooled.dn.gene.df$direction <- 'downregulated'

TableS1.df <- rbind(pooled.up.gene.df,pooled.dn.gene.df)

write_xlsx(TableS1.df[,c('gene.id','cancer.type','direction','log2FoldChange','pvalue','padj')],  "~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Table/TableS1.xlsx")


#Table S2

load('client-side/output/DE.gene.overview.R.output/DE.gene.overview.RData')

View(GO.rs.list$BRCA.Basal$up)



# Table S3
load('client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData')
driver.gene.df    <- survival.df[survival.df$p.val < 0.05 & survival.df$hr > 0 & survival.df$gene.type == 'up',]
suppresor.gene.df <- survival.df[survival.df$p.val < 0.05 & survival.df$hr < 0 & survival.df$gene.type == 'dn',]


# require(dplyr)
# library(readxl)
# library(writexl)
# 
# load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
# 
# 
# #####################################################################
# # Table S1
# ####################################################################
# 
# cut.off <- 0.05
# 
# deseq2.res         <- de.res.metastasis.liver.vs.breast.lumb
# up.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange > 1  & deseq2.res$padj < cut.off]
# dn.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange < -1 & deseq2.res$padj < cut.off]
# lumb.up            <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv", stringsAsFactors=FALSE)$x
# lumb.dn            <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.dn.csv", stringsAsFactors=FALSE)$x
# lumb.df            <- data.frame(deseq2.up= length(up.gene), deseq2.dn= length(dn.gene),pipeline.up= length(lumb.up),pipeline.dn=length(lumb.dn)    )
# 
# deseq2.res          <- de.res.metastasis.liver.vs.breast.basal
# up.gene             <- rownames(deseq2.res)[deseq2.res$log2FoldChange > 1  & deseq2.res$padj < cut.off]
# dn.gene             <- rownames(deseq2.res)[deseq2.res$log2FoldChange < -1 & deseq2.res$padj < cut.off]
# basal.up            <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
# basal.dn            <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.dn.csv", stringsAsFactors=FALSE)$x
# basal.df            <- data.frame(deseq2.up= length(up.gene), deseq2.dn= length(dn.gene),pipeline.up= length(basal.up),pipeline.dn=length(basal.dn)    )
# 
# 
# deseq2.res          <- de.res.metastasis.liver.vs.breast.her2
# up.gene             <- rownames(deseq2.res)[deseq2.res$log2FoldChange > 1  & deseq2.res$padj < cut.off]
# dn.gene             <- rownames(deseq2.res)[deseq2.res$log2FoldChange < -1 & deseq2.res$padj < cut.off]
# her2.up             <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv", stringsAsFactors=FALSE)$x
# her2.dn             <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.dn.csv", stringsAsFactors=FALSE)$x
# her2.df             <- data.frame(deseq2.up= length(up.gene), deseq2.dn= length(dn.gene),pipeline.up= length(her2.up),pipeline.dn=length(her2.dn)    )
# 
# tableS1.df          <- rbind(lumb.df,basal.df)
# tableS1.df          <- rbind(tableS1.df,her2.df)
# 
# tableS1.df$subtype   <- c('LuminalB','Basal-like','Her2-enriched')
# colnames(tableS1.df) <- c('number.of.upregulated.gene-DESeq2','number.of.downregulated.gene-DESeq2','number.of.upregulated.gene.kept.by.pipeline','number.of.downregulated.gene.kept.by.pipeline','subtype')
# tableS1.df$fraction.of.kept.DE.gene <- (tableS1.df$number.of.upregulated.gene.kept.by.pipeline + tableS1.df$number.of.downregulated.gene.kept.by.pipeline) / (tableS1.df$`number.of.upregulated.gene-DESeq2` + tableS1.df$`number.of.downregulated.gene-DESeq2`) 
# 
# write_xlsx(tableS1.df[,c('subtype','number.of.upregulated.gene-DESeq2','number.of.downregulated.gene-DESeq2','number.of.upregulated.gene.kept.by.pipeline','number.of.downregulated.gene.kept.by.pipeline','fraction.of.kept.DE.gene')],  "~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Table/TableS1.xlsx")
# mean(tableS1.df$fraction.of.kept.DE.gene)
# 
# 
# 
# #####################################################################
# # Table S3 GO enrichment analysis results
# ####################################################################
# load('client-side/output/analyze.DE.gene.R.output/analyze.DE.gene.RData')
# 
# 
# #####################################################################
# # Table S4 metastasis driver and surpressor genes
# ####################################################################
# 
# 
# load('client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData')
# 
# flag    <- basal.survival.rs.up$p.value < 0.05 & basal.survival.rs.up$effect.size > 0
# basal.survival.rs.up[flag,] %>% View
# flag    <- basal.survival.rs.dn$p.value < 0.05 & basal.survival.rs.dn$effect.size < 0
# basal.survival.rs.dn[flag,] %>% View
# 
# 
# flag    <- lumb.survival.rs.up$p.value < 0.05 & lumb.survival.rs.up$effect.size > 0
# lumb.survival.rs.up[flag,] %>% View
# flag    <- lumb.survival.rs.dn$p.value < 0.05 & lumb.survival.rs.dn$effect.size < 0
# lumb.survival.rs.dn[flag,] %>% View
# 
# 
# flag    <- her2.survival.rs.up$p.value < 0.05 & her2.survival.rs.up$effect.size > 0
# her2.survival.rs.up[flag,] %>% View
# flag    <- her2.survival.rs.dn$p.value < 0.05 & her2.survival.rs.dn$effect.size < 0
# her2.survival.rs.dn[flag,] %>% View



