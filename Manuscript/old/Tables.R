require(dplyr)
library(readxl)
library(writexl)

load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')


#####################################################################
# Table S1
####################################################################

cut.off <- 0.05

deseq2.res         <- de.res.metastasis.liver.vs.breast.lumb
up.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange > 1  & deseq2.res$padj < cut.off]
dn.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange < -1 & deseq2.res$padj < cut.off]
lumb.up            <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv", stringsAsFactors=FALSE)$x
lumb.dn            <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.dn.csv", stringsAsFactors=FALSE)$x
lumb.df            <- data.frame(deseq2.up= length(up.gene), deseq2.dn= length(dn.gene),pipeline.up= length(lumb.up),pipeline.dn=length(lumb.dn)    )

deseq2.res          <- de.res.metastasis.liver.vs.breast.basal
up.gene             <- rownames(deseq2.res)[deseq2.res$log2FoldChange > 1  & deseq2.res$padj < cut.off]
dn.gene             <- rownames(deseq2.res)[deseq2.res$log2FoldChange < -1 & deseq2.res$padj < cut.off]
basal.up            <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
basal.dn            <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.dn.csv", stringsAsFactors=FALSE)$x
basal.df            <- data.frame(deseq2.up= length(up.gene), deseq2.dn= length(dn.gene),pipeline.up= length(basal.up),pipeline.dn=length(basal.dn)    )


deseq2.res          <- de.res.metastasis.liver.vs.breast.her2
up.gene             <- rownames(deseq2.res)[deseq2.res$log2FoldChange > 1  & deseq2.res$padj < cut.off]
dn.gene             <- rownames(deseq2.res)[deseq2.res$log2FoldChange < -1 & deseq2.res$padj < cut.off]
her2.up             <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv", stringsAsFactors=FALSE)$x
her2.dn             <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.dn.csv", stringsAsFactors=FALSE)$x
her2.df             <- data.frame(deseq2.up= length(up.gene), deseq2.dn= length(dn.gene),pipeline.up= length(her2.up),pipeline.dn=length(her2.dn)    )

tableS1.df          <- rbind(lumb.df,basal.df)
tableS1.df          <- rbind(tableS1.df,her2.df)

tableS1.df$subtype   <- c('LuminalB','Basal-like','Her2-enriched')
colnames(tableS1.df) <- c('number.of.upregulated.gene-DESeq2','number.of.downregulated.gene-DESeq2','number.of.upregulated.gene.kept.by.pipeline','number.of.downregulated.gene.kept.by.pipeline','subtype')
tableS1.df$fraction.of.kept.DE.gene <- (tableS1.df$number.of.upregulated.gene.kept.by.pipeline + tableS1.df$number.of.downregulated.gene.kept.by.pipeline) / (tableS1.df$`number.of.upregulated.gene-DESeq2` + tableS1.df$`number.of.downregulated.gene-DESeq2`) 

write_xlsx(tableS1.df[,c('subtype','number.of.upregulated.gene-DESeq2','number.of.downregulated.gene-DESeq2','number.of.upregulated.gene.kept.by.pipeline','number.of.downregulated.gene.kept.by.pipeline','fraction.of.kept.DE.gene')],  "~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Table/TableS1.xlsx")
mean(tableS1.df$fraction.of.kept.DE.gene)



#####################################################################
# Table S3 GO enrichment analysis results
####################################################################
load('client-side/output/analyze.DE.gene.R.output/analyze.DE.gene.RData')


#####################################################################
# Table S4 metastasis driver and surpressor genes
####################################################################


load('client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData')

flag    <- basal.survival.rs.up$p.value < 0.05 & basal.survival.rs.up$effect.size > 0
basal.survival.rs.up[flag,] %>% View
flag    <- basal.survival.rs.dn$p.value < 0.05 & basal.survival.rs.dn$effect.size < 0
basal.survival.rs.dn[flag,] %>% View


flag    <- lumb.survival.rs.up$p.value < 0.05 & lumb.survival.rs.up$effect.size > 0
lumb.survival.rs.up[flag,] %>% View
flag    <- lumb.survival.rs.dn$p.value < 0.05 & lumb.survival.rs.dn$effect.size < 0
lumb.survival.rs.dn[flag,] %>% View


flag    <- her2.survival.rs.up$p.value < 0.05 & her2.survival.rs.up$effect.size > 0
her2.survival.rs.up[flag,] %>% View
flag    <- her2.survival.rs.dn$p.value < 0.05 & her2.survival.rs.dn$effect.size < 0
her2.survival.rs.dn[flag,] %>% View



