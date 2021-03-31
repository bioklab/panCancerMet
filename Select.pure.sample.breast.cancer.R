# Aim1: PAM50 subtyping on TCGA breast cancer samples
# Aim2: select pure TCGA breast cancer samples based on sample-cell.line correlation analysis

require(plyr)
require(dplyr)
require(genefu)
require(Rtsne)
require(ggplot2)
require(dplyr)
source('client-side/code/util.R')
load('server-side/RData/Breast Invasive Carcinoma.RData')


# TCGA.breast.cancer.polyA.sample           <- sample.meta.df$sample.id[sample.meta.df$primary.disease.or.tissue == 'Breast Invasive Carcinoma']
# TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix
# TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix

brca.tcga.meta <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/cBioPortal/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt", stringsAsFactors=FALSE,comment.char = '#',header = TRUE)
duc.sample     <- brca.tcga.meta$SAMPLE_ID[brca.tcga.meta$CANCER_TYPE_DETAILED == 'Breast Invasive Ductal Carcinoma'] %>% as.character()
duc.sample     <- intersect(duc.sample,colnames(log2.read.count.matrix))
TCGA.breast.cancer.polyA.sample           <- duc.sample
TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix[,duc.sample]
TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix[,duc.sample]


#### Perform pam50 subtyping for TCGA polyA samples###########
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') # well, I borrow some information from the metastatic breast cancer evaluation project 
pam50.gene                  <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- TCGA.breast.cancer.log2.fpkm.matrix[pam50.gene.df$probe %>% as.character,TCGA.breast.cancer.polyA.sample] %>% t 
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )


TCGA.breast.cancer.polyA.LumB.sample    <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'LumB']
TCGA.breast.cancer.polyA.Basal.sample   <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'Basal']
TCGA.breast.cancer.polyA.Her2.sample    <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'Her2']
TCGA.breast.cancer.polyA.LumA.sample    <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'LumA']
TCGA.breast.cancer.polyA.Normal.sample  <- TCGA.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'Normal']




############## TCGA.sample - CCLE.cell.line correlation analysis ################
load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000   <- names(sort(rank.sd,decreasing =TRUE))[1:1000]

#### TCGA #######
cor.cut.off                           <- 0.3
rank.cut.off                          <- 1010

rs.TCGA                               <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.Basal.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                        <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.breast.cancer.Basal.sample   <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
TCGA.best.cell.line.Basal             <- rs.TCGA$best.cell.line


rs.TCGA                               <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.Her2.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                        <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.breast.cancer.Her2.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
TCGA.best.cell.line.Her2              <- rs.TCGA$best.cell.line


rs.TCGA                               <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.LumB.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                        <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.breast.cancer.LumB.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
TCGA.best.cell.line.LumB               <- rs.TCGA$best.cell.line


rs.TCGA                               <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.LumA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                        <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.breast.cancer.LumA.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
TCGA.best.cell.line.LumA              <- rs.TCGA$best.cell.line




######## MET500 ##########

load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')

rs.MET500                             <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Basal.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                        <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.breast.cancer.Basal.sample   <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
biopsy.site                           <- MET500.sample.meta[pure.MET.breast.cancer.Basal.sample,'biopsy.site'] %>% as.character()
pure.MET.breast.cancer.Basal.sample   <- pure.MET.breast.cancer.Basal.sample[biopsy.site == 'LIVER']
MET500.best.cell.line.Basal           <- rs.MET500$best.cell.line




rs.MET500                             <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumB.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                        <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.breast.cancer.LumB.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
biopsy.site                           <- MET500.sample.meta[pure.MET.breast.cancer.LumB.sample,'biopsy.site'] %>% as.character()
pure.MET.breast.cancer.LumB.sample    <- pure.MET.breast.cancer.LumB.sample[biopsy.site == 'LIVER']
MET500.best.cell.line.LumB            <- rs.MET500$best.cell.line


rs.MET500                             <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Her2.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                        <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.breast.cancer.Her2.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
biopsy.site                           <- MET500.sample.meta[pure.MET.breast.cancer.Her2.sample,'biopsy.site'] %>% as.character()
pure.MET.breast.cancer.Her2.sample    <- pure.MET.breast.cancer.Her2.sample[biopsy.site == 'LIVER']
MET500.best.cell.line.Her2            <- rs.MET500$best.cell.line


rs.MET500                             <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                              <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                        <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.breast.cancer.LumA.sample    <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
biopsy.site                           <- MET500.sample.meta[pure.MET.breast.cancer.LumA.sample,'biopsy.site'] %>% as.character()
pure.MET.breast.cancer.LumA.sample    <- pure.MET.breast.cancer.LumA.sample[biopsy.site == 'LIVER']
MET500.best.cell.line.LumA            <- rs.MET500$best.cell.line


# attentation! SRR4238348 and SRR4305661 are duplicated samples. It is NOT our fault, maybe the author upload the wrong sequence files
pure.MET.breast.cancer.LumB.sample <- setdiff(pure.MET.breast.cancer.LumB.sample,'SRR4238348')

save(file = 'client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData',
     list = c(
              'pure.PRI.breast.cancer.LumB.sample',  'pure.PRI.breast.cancer.LumA.sample',
              'pure.PRI.breast.cancer.Her2.sample',  'pure.PRI.breast.cancer.Basal.sample',
              'pure.MET.breast.cancer.LumB.sample',  'pure.MET.breast.cancer.LumA.sample',
              'pure.MET.breast.cancer.Her2.sample',  'pure.MET.breast.cancer.Basal.sample',
              'TCGA.best.cell.line.Basal','TCGA.best.cell.line.LumB',
              'TCGA.best.cell.line.Her2' ,'TCGA.best.cell.line.LumA',
              'MET500.best.cell.line.Basal' ,'MET500.best.cell.line.LumB',
              'MET500.best.cell.line.Her2'  ,'MET500.best.cell.line.LumA'
              
              )
)





# ############### t-SNE visualization ###############
# dist.obj       <- as.dist(1- cor(log2.fpkm.matrix[pam50.gene,TCGA.breast.cancer.polyA.sample],method='spearman'))
# 
# set.seed(8) # I want to reproduce the tsne results, 8 is just a arbitrary numnber, it DOES NOT change the conclusion
# tsne.rs        <- Rtsne(dist.obj,perplexity = 15)
# TCGA.pam50.subtype.vec <- c(   rep(x='LuminalA',       times= TCGA.breast.cancer.polyA.LumA.sample %>% length),
#                                rep(x='LuminalB',       times= TCGA.breast.cancer.polyA.LumB.sample %>% length),
#                                rep(x='Her2-enriched',  times= TCGA.breast.cancer.polyA.Her2.sample  %>% length),
#                                rep(x='Basal-like',     times= TCGA.breast.cancer.polyA.Basal.sample %>% length),
#                                rep(x='Normal-like',    times= TCGA.breast.cancer.polyA.Normal.sample %>% length)
#                            )
# names(TCGA.pam50.subtype.vec) <- c(  TCGA.breast.cancer.polyA.LumA.sample,
#                                      TCGA.breast.cancer.polyA.LumB.sample,
#                                      TCGA.breast.cancer.polyA.Her2.sample,
#                                      TCGA.breast.cancer.polyA.Basal.sample,
#                                      TCGA.breast.cancer.polyA.Normal.sample
#                                    )
# TCGA.pam50.subtype.tsne.df <- data.frame(dim1=tsne.rs$Y[,1],
#                       dim2=tsne.rs$Y[,2],
#                       subtype       = TCGA.pam50.subtype.vec[TCGA.breast.cancer.polyA.sample]
# )
# ggplot(TCGA.pam50.subtype.tsne.df,aes(x=dim1,y=dim2,color=subtype)) + geom_point(size=6) +  scale_shape_manual(values=c(8,15:18))
# 



# save(file = 'client-side/output/TCGA.breast.cancer.meta.R.output/TCGA.breast.cancer.meta.RData',
#      list = c('TCGA.breast.cancer.polyA.LumB.sample',  'TCGA.breast.cancer.polyA.LumA.sample',
#               'TCGA.breast.cancer.polyA.Her2.sample',  'TCGA.breast.cancer.polyA.Basal.sample',
#               'TCGA.breast.cancer.polyA.Normal.sample',
#               'pure.TCGA.breast.cancer.polyA.LumB.sample',  'pure.TCGA.breast.cancer.polyA.LumA.sample',
#               'pure.TCGA.breast.cancer.polyA.Her2.sample',  'pure.TCGA.breast.cancer.polyA.Basal.sample',
#               'TCGA.pam50.subtype.tsne.df'
#      )
# )










############### Trash code #################
# tt.Basal <- cancer.data[,TCGA.breast.cancer.polyA.Basal.sample]
# tt.LumA  <- cancer.data[,TCGA.breast.cancer.polyA.LumA.sample]
# tt.LumB  <- cancer.data[,TCGA.breast.cancer.polyA.LumB.sample]
# tt.Her2 <- cancer.data[,TCGA.breast.cancer.polyA.Her2.sample]
# 
# g.vec <- c('ENSG00000186847','ENSG00000111057','ENSG00000124107','ENSG00000148513')
# 
# data <- cbind(tt.LumA,tt.LumB,tt.Basal,tt.Her2)
# pca.rs <- prcomp(data[g.vec,] %>% t)
# plot(pca.rs$x[,1:2])
# 
# tsne.rs <- Rtsne(data[g.vec,] %>% t)
# 
# draw.df <- data.frame(pca.rs$x,tsne.rs$Y)
# draw.df$subtype <- 'Basal'
# draw.df[rownames(draw.df) %in% TCGA.breast.cancer.polyA.LumA.sample,'subtype'] <- 'LumA'
# draw.df[rownames(draw.df) %in% TCGA.breast.cancer.polyA.LumB.sample,'subtype'] <- 'LumB'
# draw.df[rownames(draw.df) %in% TCGA.breast.cancer.polyA.Her2.sample,'subtype'] <- 'Her2'
# 
# ggplot(draw.df)+geom_point(aes(x=X1,y=X2,color=subtype))

# #Well, let us run tSNE with the 1000 most varied genes
# g <- intersect(CCLE.rna.seq.marker.gene.1000,rownames(log2.fpkm.matrix))
# 
#load('~/Project/InSilicoCRISPR/client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData')
# #dist.obj       <- as.dist(1- cor(log2.fpkm.matrix[g,TCGA.breast.cancer.polyA.sample],method='spearman'))
# purity         <- TCGA.CD8.T.cell.level.df[TCGA.breast.cancer.polyA.sample,'tumor.purity']
# names(purity)  <- TCGA.breast.cancer.polyA.sample
# h.sample <- names(purity)[purity >= 0.75]
# h.sample <- h.sample[is.na(h.sample) == FALSE]
# h.sample <- setdiff(h.sample,TCGA.breast.cancer.polyA.Basal.sample)
# g <- intersect(CCLE.rna.seq.marker.gene.1000,rownames(log2.fpkm.matrix))
# 
# dist.obj       <- as.dist(1- cor(log2.fpkm.matrix[g,h.sample],method='spearman'))
# set.seed(8) # I want to reproduce the tsne results, 8 is just a arbitrary numnber, it DOES NOT change the conclusion
# tsne.rs        <- Rtsne(dist.obj,perplexity = 8)
# pca.rs <- prcomp(log2.fpkm.matrix[g,h.sample] %>% t)
# 
# draw.df <- data.frame(tsne1=tsne.rs$Y[,1],
#                       tsne2=tsne.rs$Y[,2],
#                       pc1=pca.rs$x[,1],
#                       pc2=pca.rs$x[,6],
#                       subtype       = TCGA.pam50.subtype.vec[h.sample]
# )
# tt.sampe <- rownames(draw.df)[draw.df$tsne2 < -25]
# uu.sample <- setdiff(h.sample,tt.sample)
# draw.df$is.APLP1 <- ifelse(rownames(draw.df) %in% tt.sample,'yes','no')
# ggplot(draw.df,aes(x=tsne1,y=tsne2,color=is.APLP1)) + geom_point(size=6) +  scale_shape_manual(values=c(8,15:18))
# 
# rs.TCGA       <- pick.out.cell.line(expr.of.samples = log2.fpkm.matrix[,tt.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# 
# pairs(pca.rs$x[,c(1,6)],col=ifelse(rownames(pca.rs$x) %in% tt.sample,'red','black'))
# 
# pca.rs <- prcomp(log2.fpkm.matrix[g,c(tt.sample,h.sample[17:33])] %>% t)
# 
# #APLP1
# APLP1 <- 'ENSG00000105290'
# uu.sample <- setdiff(h.sample,tt.sample)
# flag <- apply(log2.fpkm.matrix[,h.sample],1,function(x) median(x) > log2(1+0.1))
# e.gene <- rownames(log2.fpkm.matrix)[flag]
# 
# tt.matrix <- log2.fpkm.matrix[e.gene,tt.sample]
# uu.matrix <- log2.fpkm.matrix[e.gene,uu.sample]
# delta     <- apply(tt.matrix,1,median) -  apply(uu.matrix,1,median)
# p.value   <- sapply(1:nrow(tt.matrix),function(i) wilcox.test(tt.matrix[i,],uu.matrix[i,])$p.value )
# rs        <- data.frame(effect.size=delta,p.value=p.value)
# rownames(rs) <- e.gene
# rs$fdr <- p.adjust(rs$p.value,method='fdr')
# 
# up.gene <- rownames(rs)[rs$fdr < 0.01 & rs$effect.size > 1.5]
# dn.gene <- rownames(rs)[rs$fdr < 0.01 & rs$effect.size < -1.5]
# write.csv(x=up.gene,file='~/Desktop/up.gene.csv',quote=FALSE)
# write.csv(x=dn.gene,file='~/Desktop/dn.gene.csv',quote=FALSE)
# 
# up.gene.rs <- (rs)[rs$fdr < 0.01 & rs$effect.size > 1,]
# dn.gene.rs <- (rs)[rs$fdr < 0.01 & rs$effect.size < -1,]
# up.gene.rs <- up.gene.rs[order(up.gene.rs$effect.size,decreasing = TRUE),]
# dn.gene.rs <- dn.gene.rs[order(dn.gene.rs$effect.size,decreasing = TRUE),]
# 
# rs.tt.sample     <- pick.out.cell.line(expr.of.samples = log2.fpkm.matrix[,tt.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# 
# 
# rs <- foreach(g=e.gene,.combine='rbind') %do% {
#    delta <- log2.fpkm.matrix[g,tt.sample] %>% median  -  log2.fpkm.matrix[g,uu.sample] %>% median
#    p.value <- wilcox.test(log2.fpkm.matrix[g,tt.sample],log2.fpkm.matrix[g,uu.sample])$p.value  
#    data.frame(effect.size=delta,p.value=p.value)
# }
# rownames(rs) <- e.gene

# get.TCGA.cohort.id <- function(x) {
#   tmp <- strsplit(x,split = '\\-') %>% unlist  
#   paste(tmp[1:3],collapse = '-')
#   
# }
# 
# 
# TCGA.breast.cancer.polyA.sample         <- colnames(log2.fpkm.matrix)
# TCGA.breast.cancer.polyA.cohort         <- sapply(TCGA.breast.cancer.polyA.sample,get.TCGA.cohort.id)
# tmp                                     <- TCGA.clinical.data.list$BRCA
# Ductal.carcinoma.cohort                 <- rownames(tmp)[tmp$histological_type == 'infiltrating ductal carcinoma']
# TCGA.breast.cancer.polyA.sample         <- TCGA.breast.cancer.polyA.sample[TCGA.breast.cancer.polyA.cohort %in% Ductal.carcinoma.cohort]


