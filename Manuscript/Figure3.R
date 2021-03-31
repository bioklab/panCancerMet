require(ggplot2)
require(CePa)
require(dplyr)
require(ComplexHeatmap)
require(RColorBrewer)
require(circlize)
require(segmented)
source('client-side/code/Manuscript/ggplot.style.R')
require(foreach)

################################################################################################################################################################      
#----------------------------- Figure 3 ---------------------------------------#
################################################################################################################################################################     

##############################################################################
# (a): GO enrichment analysis results
##############################################################################

load('client-side/output/DE.gene.overview.R.output/DE.gene.overview.RData')
most.sig.GO.rs.list <- foreach(GO.rs = GO.rs.list) %do% {
  list(up = GO.rs$up[1,], dn = GO.rs$dn[1,])
}
names(most.sig.GO.rs.list) <- names(GO.rs.list)

df <- foreach(item = most.sig.GO.rs.list,.combine='rbind') %do% {
  item$up
}
df$cancer.type <- names(most.sig.GO.rs.list)
up.df          <- df
up.df$class    <- 'up'

df <- foreach(item = most.sig.GO.rs.list,.combine='rbind') %do% {
  item$dn
}
df$cancer.type <- names(most.sig.GO.rs.list)
dn.df          <- df
dn.df$class    <- 'dn'

draw.df <- rbind(up.df,dn.df)
draw.df$cancer.type <- factor(draw.df$cancer.type, levels = c('BRCA.Basal','BRCA.LumB','BRCA.Her2','PRAD','COAD','PNET','SINET'))
ggplot(draw.df,aes(x=cancer.type,y= -1 * log10(pvalue),fill=class))  + geom_bar(stat="identity",position="dodge",width= 0.7,show.legend = FALSE) + 
  coord_flip() + theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                                  axis.text  = element_text( size=25, face="bold"),
                                                  plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                                  axis.line.x = element_line(colour = "black",size = 3),
                                                  axis.line.y = element_line(colour = "black",size = 3))  + 
  scale_fill_manual(values = c('up' = 'red', 'dn' = 'blue'))
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/GO.bar.plot.pdf',width = 20,height=20)




################################################################################   
# (b): boxplot of G1/S ssGSEA scores
################################################################################   
load('client-side/output/cell.cycle.R.output/cell.cycle.RData')

ggplot(cell.cycle.score.df[cell.cycle.score.df$phase == 'G1/S',], aes(x=cancer.type, y=score, fill=site)) + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(),size=3.0,show.legend = FALSE) + 
  ylab('G1/S ssGSEA score') + ggplot.style  + scale_fill_manual(values = c('primary'='#EF8A62','metastatic'='#D1E5F0')) + ylim(-7,7)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/G1.and.S.score.pdf',width = 20,height=20)

test.diff <- function(phase,cancer.type) {
    f1 <-   cell.cycle.score.df$phase == phase & cell.cycle.score.df$cancer.type == cancer.type & cell.cycle.score.df$site == 'primary'
    f2 <-   cell.cycle.score.df$phase == phase & cell.cycle.score.df$cancer.type == cancer.type & cell.cycle.score.df$site == 'metastatic'
    wilcox.test(cell.cycle.score.df$score[f1],cell.cycle.score.df$score[f2])$p.value
}

G1.and.S.pval <- sapply(cell.cycle.score.df$cancer.type %>% unique %>% as.character, function(x) test.diff('G1/S',x))
G1.and.S.fdr  <- p.adjust(G1.and.S.pval,method='fdr')
G2.and.M.pval <- sapply(cell.cycle.score.df$cancer.type %>% unique %>% as.character, function(x) test.diff('G2/M',x))
G2.and.M.fdr  <- p.adjust(G2.and.M.pval,method='fdr')


################################################################################   
# (c) : boxplot of G2/M - G1/S ssGSEA score
################################################################################  

ggplot(phase.score.diff.df,aes(x=cancer.type,y=score,fill=site)) + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(),size=3.0,show.legend = FALSE) + 
  ylab('ssGSEA score difference(G2/M - G1/S)') + ggplot.style  + scale_fill_manual(values = c('primary'='#EF8A62','metastatic'='#D1E5F0')) + ylim(-7,7)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/phase.diff.score.pdf',width = 20,height=20)


test.diff <- function(cancer.type) {
  f1 <-    phase.score.diff.df$cancer.type == cancer.type & phase.score.diff.df$site == 'primary'
  f2 <-    phase.score.diff.df$cancer.type == cancer.type & phase.score.diff.df$site == 'metastatic'
  wilcox.test(phase.score.diff.df$score[f1],phase.score.diff.df$score[f2])$p.value
}
p.value <- sapply(phase.score.diff.df$cancer.type %>% unique %>% as.character, function(x) test.diff(x))











################################################################################################################################################################      
#----------------------------- Figure Sx ---------------------------------------#
################################################################################################################################################################     

################################################################################   
# (a): 
################################################################################  

ggplot(cell.cycle.score.df[cell.cycle.score.df$phase == 'G2/M',], aes(x=cancer.type, y=score, fill=site)) + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(),size=3.0,show.legend = FALSE) + 
  ylab('G2/M ssGSEA score') + ggplot.style  + scale_fill_manual(values = c('primary'='#EF8A62','metastatic'='#D1E5F0')) + ylim(-7,7)

ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/G2.and.M.score.pdf',width = 20,height=20)



################################################################################   
# (b): VennDiagram and GO enrichment 
################################################################################  

load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')
load('client-side/output/DE.prostate.cancer.SRP253428.R.output/DE.prostate.cancer.SRP253428.RData')
require(clusterProfiler)
library(org.Hs.eg.db)


gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id


perform.GO.analysis <- function(DE.rs) {
  up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = intersect(DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,mapping.df$ensembl_gene_id))
  up.gene.annotation.df  <- up.gene.annotation.df[complete.cases(up.gene.annotation.df),]
  GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL %>% unique ,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
  GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 , ]
  GO.rs.1.up             <- GO.rs.1.up[order(GO.rs.1.up$pvalue),]
  
  dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = intersect(DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene,mapping.df$ensembl_gene_id))
  dn.gene.annotation.df  <- dn.gene.annotation.df[complete.cases(dn.gene.annotation.df),]
  GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL %>% unique,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
  GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 , ]
  GO.rs.1.dn             <- GO.rs.1.dn[order(GO.rs.1.dn$pvalue),]
  
  list(up = GO.rs.1.up,dn = GO.rs.1.dn)
}

PRAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene %>% length
PRAD.SRP253428.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene %>% length


ov.up.gene        <- intersect(PRAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,PRAD.SRP253428.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene) 
length(ov.up.gene)

ov.up.gene.symbol <- mapping.df[ov.up.gene,'symbol'] %>% unique
GO.ov.up.gene     <- enrichGO(gene=ov.up.gene.symbol ,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result


################################################################################   
# (c,d,e):  valadition with SRP253428 (aka GSE147250) dataset (for PRAD)
################################################################################  

load('client-side/output/cell.cycle.R.output/SRP253428.validation.RData')

pri.flag <- G1.and.S.scroe.df$site == 'primary'
met.flag <- G1.and.S.scroe.df$site == 'metastatic'


ggplot(G1.and.S.scroe.df, aes(x=site,y= scale(score), fill=site)) + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(),size=3.0,show.legend = FALSE) + 
  ylab('G1/S ssGSEA score') + ggplot.style  + scale_fill_manual(values = c('primary'='#EF8A62','metastatic'='#D1E5F0')) + ylim(-7,7)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/SRP253428.G1.and.S.score.pdf',width = 20,height=20)
wilcox.test(G1.and.S.scroe.df$score[pri.flag], G1.and.S.scroe.df$score[met.flag])


ggplot(G2.and.M.scroe.df, aes(x=site,y= scale(score), fill=site)) + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(),size=3.0,show.legend = FALSE) + 
  ylab('G2/M ssGSEA score') + ggplot.style  + scale_fill_manual(values = c('primary'='#EF8A62','metastatic'='#D1E5F0')) + ylim(-7,7)
ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/SRP253428.G2.and.M.score.pdf',width = 20,height=20)
wilcox.test(G2.and.M.scroe.df$score[pri.flag], G2.and.M.scroe.df$score[met.flag])



ggplot(phase.diff.score.df, aes(x=site,y= scale(score), fill=site)) + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(),size=3.0,show.legend = FALSE) + 
  ylab('G2/M - G1/S ssGSEA score') + ggplot.style  + scale_fill_manual(values = c('primary'='#EF8A62','metastatic'='#D1E5F0')) + ylim(-7,7)

ggsave(filename = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/SRP253428.phase.diff.score.pdf',width = 20,height=20)
wilcox.test(phase.diff.score.df$score[pri.flag], phase.diff.score.df$score[met.flag])


# ##############################################
# # Fig 3a: expressin pattern of upregulated tumor-intrinsic DE genes (LuminalB)
# ##############################################
# median.tpm.matrix <- read.gct(file='client-side/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
# get.gene.id <- function(x) {
#   strsplit(x=x,split='\\.') %>% unlist %>% head(1)
# }
# rownames(median.tpm.matrix) <- sapply(rownames(median.tpm.matrix),get.gene.id)
# 
# 
# draw.color.bar <- function(value.range, color.range){
#   f <- colorRamp2(value.range,color.range)
#   x <- sapply(value.range, f)
#   plot(c(0, 40*length(x)), c(0, 100), type = "n")
#   i <- 0:(length(x)-1)
#   rect(i*40, 0, (i+1)*40, 100, col=x)
# }
# pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/color.bar.pdf',width = 20,height=20)
# draw.color.bar(value.range = seq(from=-5,to=5,by=0.5),
#                color.range = colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev
# )
# dev.off()
# 
# lumb.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv",  stringsAsFactors=FALSE)$x
# her2.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv",  stringsAsFactors=FALSE)$x
# basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
# 
# lumb.gene  <- intersect(c(lumb.up.gene),  rownames(median.tpm.matrix))
# her2.gene  <- intersect(c(her2.up.gene),  rownames(median.tpm.matrix))
# basal.gene <- intersect(c(basal.up.gene), rownames(median.tpm.matrix))
# 
# 
# tmp           <- apply((median.tpm.matrix[lumb.gene,]+1) %>% log2,1,scale) %>% t
# colnames(tmp) <- colnames(median.tpm.matrix)
# flag          <- which(rowSums(tmp) > 0)
# pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/lumb.tumor.intrinsic.upregulated.DE.gene.heatmap.pdf',width = 20,height=20)
# Heatmap(tmp[flag,], col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = TRUE,show_heatmap_legend = FALSE,cluster_columns =  FALSE)
# dev.off()
# 
# 
# 
# 
# ##############################################
# # Fig 3b: comparision GO annotation results before/after DEBoost
# ##############################################
# load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
# 
# 
# perform.GO.enrichment <- function(DEBoost.gene,DESeq2.gene){
#     DEBoost.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = DEBoost.gene)
#     DEBoost.GO.BP               <- enrichGO(gene=DEBoost.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
#     DEBoost.GO.BP               <- DEBoost.GO.BP[DEBoost.GO.BP$p.adjust < 0.05 & DEBoost.GO.BP$Count >=5,]             
#   
#     DESeq2.gene.annotation.df   <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = DESeq2.gene)
#     DESeq2.GO.BP                <- enrichGO(gene=DESeq2.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
#     DESeq2.GO.BP                <- DESeq2.GO.BP[DESeq2.GO.BP$p.adjust < 0.05 & DESeq2.GO.BP$Count >=5,]             
# 
#     GO.gene.symbol <- sapply(DESeq2.GO.BP$geneID , function(x) strsplit(x,split='/') %>% unlist) %>% unlist %>% unique
#     r              <- 1- intersect(GO.gene.symbol,DEBoost.gene.annotation.df$SYMBOL) %>% length / length(GO.gene.symbol)
#     list(DEBoost.GO.BP=DEBoost.GO.BP,DESeq2.GO.BP=DESeq2.GO.BP,uDE.gene.ratio = r)
# }
#   
# # Basal-like, upregulated gene
# DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
# flag                        <- de.res.metastasis.liver.vs.breast.basal$log2FoldChange > 1 & de.res.metastasis.liver.vs.breast.basal$padj < 0.05
# DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.basal)[flag]
# basal.up.rs                 <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)
# 
# # Basal-like, dnregulated gene
# DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.dn.csv", stringsAsFactors=FALSE)$x
# flag                        <- de.res.metastasis.liver.vs.breast.basal$log2FoldChange < -1 & de.res.metastasis.liver.vs.breast.basal$padj < 0.05
# DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.basal)[flag]
# basal.dn.rs                 <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)
# 
# 
# # Lumb, upregulated gene
# DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv", stringsAsFactors=FALSE)$x
# flag                        <- de.res.metastasis.liver.vs.breast.lumb$log2FoldChange > 1 & de.res.metastasis.liver.vs.breast.lumb$padj < 0.05
# DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.lumb)[flag]
# lumb.up.rs                  <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)
# 
# # Lumb, dnregulated gene
# DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.dn.csv", stringsAsFactors=FALSE)$x
# flag                        <- de.res.metastasis.liver.vs.breast.lumb$log2FoldChange < -1 & de.res.metastasis.liver.vs.breast.lumb$padj < 0.05
# DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.lumb)[flag]
# lumb.dn.rs                  <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)
# 
# 
# # Her2, upregulated gene
# DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv", stringsAsFactors=FALSE)$x
# flag                        <- de.res.metastasis.liver.vs.breast.her2$log2FoldChange > 1 & de.res.metastasis.liver.vs.breast.her2$padj < 0.05
# DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.her2)[flag]
# her2.up.rs                  <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)
# 
# # Her2, dnregulated gene
# DEBoost.gene                <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.dn.csv", stringsAsFactors=FALSE)$x
# flag                        <- de.res.metastasis.liver.vs.breast.her2$log2FoldChange < -1 & de.res.metastasis.liver.vs.breast.her2$padj < 0.05
# DESeq2.gene                 <- rownames(de.res.metastasis.liver.vs.breast.her2)[flag]
# her2.dn.rs                  <- perform.GO.enrichment(DEBoost.gene = DEBoost.gene,DESeq2.gene = DESeq2.gene)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ################################################################################   
# #----------------------------- Figure S3 ---------------------------------------#
# ################################################################################  
# 
# ##############################################
# # Fig S3a, S3b: expressin pattern of upregulated tumor-intrinsic DE genes (LuminalB)
# ##############################################
# tmp           <- apply((median.tpm.matrix[her2.gene,]+1) %>% log2,1,scale) %>% t
# colnames(tmp) <- colnames(median.tpm.matrix)
# flag          <- which(rowSums(tmp) > 0)
# pdf(file= '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/her2.tumor.intrinsic.upregulated.DE.gene.heatmap.pdf',width = 20,height=20)
# Heatmap(tmp[flag,], col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = TRUE,show_heatmap_legend = FALSE,cluster_columns =  FALSE)
# dev.off()
# 
# tmp           <- apply((median.tpm.matrix[basal.gene,]+1) %>% log2,1,scale) %>% t
# colnames(tmp) <- colnames(median.tpm.matrix)
# flag          <- which(rowSums(tmp) > 0)
# pdf(file = '~/OneDrive - Michigan State University/Project/BreastCancerMetaPotenial/Manuscript/Manuscript/section.Figure/Section3/basal.tumor.intrinsic.upregulated.DE.gene.heatmap.pdf',width = 20,height=20)
# Heatmap(tmp[flag,], col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),show_row_names = FALSE,show_column_names = TRUE,show_heatmap_legend = FALSE,cluster_columns =  FALSE)
# dev.off()
# 


# ################################################################################   
# #Figure 3e
# ################################################################################  
# 
# 
# ggplot(GO.ov.up.gene[1:5,],aes(x=ID,y= -1 * log10(pvalue)))  + geom_bar(stat="identity",position="dodge",width= 0.7,show.legend = FALSE) + 
#   coord_flip() + theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
#                                                   axis.text  = element_text( size=25, face="bold"),
#                                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
#                                                   axis.line.x = element_line(colour = "black",size = 3),
#                                                   axis.line.y = element_line(colour = "black",size = 3))  + 
#   scale_fill_manual(values = c('up' = 'red', 'dn' = 'blue'))
# 
# 
# 










