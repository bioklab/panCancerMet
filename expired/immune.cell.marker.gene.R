# well, maybe we should add adipocyte, marker gene CRLF1
require(plyr)
require(dplyr)

immune.cell.signature <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/immune.cell.signature.csv", comment.char="#", stringsAsFactors=FALSE)
immune.cell.signature <- immune.cell.signature[,1:2]


require(org.Hs.eg.db)
mapping.df <- AnnotationDbi::select(x=org.Hs.eg.db,keys=immune.cell.signature$gene.symbol %>% as.character,keytype = 'SYMBOL',columns = "ENSEMBL")
colnames(mapping.df) <- c('gene.symbol','ensemble.id')
immune.cell.signature <- merge(immune.cell.signature,mapping.df,all.x=TRUE)

immune.cell.signature <- immune.cell.signature[immune.cell.signature$gene.symbol %in% c('BLK','IRF7','LILRA4','PADI4') == FALSE,]

### The four genes has multiple mappings, so I first remove them and mannunaly add them back
df <- data.frame(gene.symbol= c('BLK','IRF7','LILRA4','PADI4'),
                 cell.type  = c('B cells','pDCs','pDCs','Neutrophils'),
                 ensemble.id = c('ENSG00000136573','ENSG00000185507','ENSG00000239961','ENSG00000159339')
                 )
immune.cell.signature <- rbind(immune.cell.signature,df)
immune.cell.signature[immune.cell.signature$gene.symbol == 'C15orf53','ensemble.id'] <- 'ENSG00000175779'
immune.cell.signature <- immune.cell.signature[order(immune.cell.signature$cell.type),]


immune.gene <- read.csv(file='client-side/output/analyze.DE.gene.R.output/immune.gene.list.csv')$x %>% as.character()
flag        <- immune.cell.signature$ensemble.id %in% immune.gene
immune.cell.signature <- immune.cell.signature[flag,]
save(list = c('immune.cell.signature'),file ='client-side/output/immune.cell.marker.gene.R.output/immune.cell.marker.gene.RData' )


# load('server-side//RData//Liver.RData')
# female.sample                     <- sample.meta.df$sample.id[sample.meta.df$gender == 'Female'] %>% as.character()
# Ref.liver.log2.fpkm.matrix        <- log2.fpkm.matrix[,female.sample]
# liver.median.expr                 <- apply(Ref.liver.log2.fpkm.matrix[immune.cell.signature$ensemble.id %>% as.character(),],1,median)
# 
# 
# 
# load('server-side/RData//Breast Invasive Carcinoma.RData')
# load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
# TCGA.breast.cancer.log2.fpkm.matrix <- log2.fpkm.matrix
# basal.median.expr                   <- apply(TCGA.breast.cancer.log2.fpkm.matrix[immune.cell.signature$ensemble.id %>% as.character(),pure.TCGA.breast.cancer.polyA.Basal.sample],1,median)
# lumb.median.expr                    <- apply(TCGA.breast.cancer.log2.fpkm.matrix[immune.cell.signature$ensemble.id %>% as.character(),pure.TCGA.breast.cancer.polyA.LumB.sample],1,median)
# 
# 
# load('~/Project/Cancer2CellLine/server-side/RData/CCLE_BRACA.RData')
# basal.expr <- BRACA.log2.fpkm.matrix[immune.cell.signature$ensemble.id %>% as.character(),'HCC1569']
# 
# require(xCell)
# 
# 
# xCell.corrected.gene.set <- foreach(x=xCell.data$signatures) %do% {
#   x@geneIds
# }
# names(xCell.corrected.gene.set) <- names(xCell.data$signatures) ############# remove DE genes from xCell signatures. Maybe here we should also remove genes DE between metastatic and primary
# 
# #L <- list()
# tmp.df <- foreach(i = names(xCell.corrected.gene.set),.combine='rbind') %do% {
#   tmp <- strsplit(x = i,split='%')[[1]] %>% unlist
#   cell.type  <- tmp[1]
#   data.source <- tmp[2]
#   data.frame(cell.type=cell.type,data.source=data.source,gene.id=xCell.corrected.gene.set[[i]])
# 
#   #L[[cell.type]] <- c(L[[cell.type]],xCell.corrected.gene.set[[i]])
#   #L[[cell.type]] <- unique(L[[cell.type]])
# }
# #xCell.corrected.gene.set <- L
# get.overlapped.gene <- function(x) {
#   l <- x$data.source %>% unique %>% length
#   ifelse(l >=2,1,0)
# }
# 
# gene.df <- ddply(tmp.df,.(cell.type,gene.id),get.overlapped.gene)
# gene.df <- gene.df[gene.df$V1 == 1,]
# 
# 
# 
# 
# mapping.df <- AnnotationDbi::select(x=org.Hs.eg.db,keys=gene.df$gene.id %>% as.character,keytype = 'SYMBOL',columns = "ENSEMBL")
# colnames(mapping.df) <- c('gene.symbol','ensemble.id')



# #### immune cell signature used in liver-metastasis vs breast primary comparision. Remove DE genes between liver and breast
# load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
# cut.off <- 0.05
# res     <- de.res.liver.vs.breast
# de.res.liver.vs.breast.up.gene <- rownames(res)[res$log2FoldChange > 1  & res$padj < cut.off]
# de.res.liver.vs.breast.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < cut.off]
# pooled.de.gene                 <- c(de.res.liver.vs.breast.up.gene,de.res.liver.vs.breast.dn.gene) %>% unique
# liver.vs.breast.immune.cell.signature <- immune.cell.signature[immune.cell.signature$ensemble.id %in% pooled.de.gene == FALSE,]
# 
# 
# ############## create signatures for different subtype comparison, need to remove DE genes between Basal cells and other cell types ##############
# load('server-side/RData/GSE113197_FPKM.RData')
# GSE113197 <- as.matrix(GSE113197)
# 
# require(GEOquery)
# GSE113197.obj <- getGEO('GSE113197',GSEMatrix=FALSE)
# gsm.list      <- GSE113197.obj@gsms
# 
# gsm.title.vec <- sapply(gsm.list,function(x) x@header$title)
# parse.gsm.title <- function(x) {
#     tmp <- strsplit(x = x,split = "_") %>% unlist  
#     data.frame( batch.id= paste(tmp[1:2],collapse="_"), cell.type=tmp[3])
# }
# GSE113197.meta.df <- foreach(x=gsm.title.vec,.combine = 'rbind') %do% {
#     parse.gsm.title(x)
# }
# rownames(GSE113197.meta.df) <- names(gsm.title.vec)
# 
# 
# 
# sc.expr.matrix <- log2(GSE113197 + 1)
# 
# basal.sample <- rownames(GSE113197.meta.df)[GSE113197.meta.df$batch.id == 'L1_I1' & GSE113197.meta.df$cell.type == 'BAS']
# lum.sample   <- rownames(GSE113197.meta.df)[GSE113197.meta.df$batch.id == 'L1_I1' & GSE113197.meta.df$cell.type == 'LUM']
# tmp          <- sc.expr.matrix[,c(basal.sample,lum.sample)]
# 
# 
# 
# basal.sample <- rownames(GSE113197.meta.df)[GSE113197.meta.df$batch.id == 'L1_I1' & GSE113197.meta.df$cell.type == 'BAS']
# lum.sample   <- rownames(GSE113197.meta.df)[GSE113197.meta.df$batch.id == 'L1_I1' & GSE113197.meta.df$cell.type == 'LUM']
# tmp          <- sc.expr.matrix[,c(basal.sample,lum.sample)]
# 
# 
# tsne.rs <- Rtsne::Rtsne(X=tmp %>% t)
# draw.df <- data.frame(tsne.rs$Y[,1:2],cell.type = GSE113197.meta.df[colnames(tmp),'cell.type'])
# ggplot(draw.df) + geom_point(aes(x=X1,y=X2,color=cell.type))
# 
# 
# 
# 
# g.vec <- intersect(immune.cell.signature$ensemble.id %>% as.character,rownames(GSE113197))
# p.value.vec <- foreach(g=g.vec,.combine='c') %do% {
#     wilcox.test( tmp[g,basal.sample],tmp[g,lum.sample])$p.value 
# }
# names(p.value.vec) <- g.vec
# p.value.vec        <- p.value.vec[is.na(p.value.vec) == FALSE]
# fdr.vec            <- p.adjust(p.value.vec,method='fdr')
# de.gene            <- names(fdr.vec)[fdr.vec < 0.05]
# breast.basal.subtype.immune.cell.signature <- immune.cell.signature[immune.cell.signature$ensemble.id %in% de.gene == FALSE,]
# 
# #########################
# save(list = c('liver.vs.breast.immune.cell.signature','breast.basal.subtype.immune.cell.signature','immune.cell.signature'),file ='client-side/output/immune.cell.marker.gene.R.output/immune.cell.marker.gene.RData' )
# 






################## Trash code ##############

#####################################################################
# marker.gene <- c('ENSG00000167286', #CD3D, T cell
#                  'ENSG00000010610', #CD4,  T cell
#                  'ENSG00000153563', #CD8A, T cell
#                  'ENSG00000172116', #CD8B, T cell   
#                  'ENSG00000177455', #CD19, B cell
#                  'ENSG00000156738', #CD20, B cell
#                  'ENSG00000149294', #CD56, NK cell
#                  'ENSG00000005961', #CD41, Platlet
#                  'ENSG00000140678', #CD11c,Dendritic cell
#                  'ENSG00000170458'  #CD14, Macrophage
# )
# marker.gene.symbol <- c('CD3D', #CD3D, T cell
#                         'CD4',  #CD4,  T cell
#                         'CD8A', #CD8A, T cell
#                         'CD8B', #CD8B, T cell   
#                         'CD19', #CD19, B cell
#                         'CD20', #CD20, B cell
#                         'CD56', #CD56, NK cell
#                         'CD41', #CD41, Platlet
#                         'CD11c',#CD11c,Dendritic cell
#                         'CD14'  #CD14, Macrophage
# )


# xCell.corrected.gene.set <- foreach(x=xCell.data$signatures) %do% {
#     x@geneIds
# }
# names(xCell.corrected.gene.set) <- names(xCell.data$signatures) 
# 
# tmp.df <- foreach(i = names(xCell.corrected.gene.set),.combine='rbind') %do% {
#     tmp <- strsplit(x = i,split='%')[[1]] %>% unlist      
#     cell.type  <- tmp[1]
#     data.source <- tmp[2]
#     data.frame(cell.type=cell.type,data.source=data.source,gene.id=xCell.corrected.gene.set[[i]])
# }
# 
# get.overlapped.gene <- function(x) {
#     l <- x$data.source %>% unique %>% length
#     ifelse(l >=2,1,0)
# }
# 
# gene.df        <- ddply(tmp.df,.(cell.type,gene.id),get.overlapped.gene)
# gene.df        <- gene.df[gene.df$V1 == 1,]
# immune.gene.df <- gene.df
# colnames(immune.gene.df) <- c('cell.type','gene.symbol','flag')
# 
# 
# require(org.Hs.eg.db)
# keytypes(org.Hs.eg.db)
# require(AnnotationDbi)
# mapping.df <- AnnotationDbi::select(x=org.Hs.eg.db,keys=immune.gene.df$gene.symbol %>% as.character,keytype = 'SYMBOL',columns = "ENSEMBL")
# colnames(mapping.df) <- c('gene.symbol','ensemble.id')
# immune.gene.df <- merge(immune.gene.df,mapping.df,all.x=TRUE)
# 
# 
# 
# ############# remove DE genes from xCell signatures. Maybe here we should also remove genes DE between metastatic and primary 
# lumb.up.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.R.output/lumb.up.csv", stringsAsFactors=FALSE)$x %>% as.character()
# lumb.dn.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.R.output/lumb.dn.csv", stringsAsFactors=FALSE)$x %>% as.character()
# luma.up.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.R.output/luma.up.csv", stringsAsFactors=FALSE)$x %>% as.character()
# luma.dn.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.R.output/luma.dn.csv", stringsAsFactors=FALSE)$x %>% as.character()
# her2.up.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.R.output/her2.up.csv", stringsAsFactors=FALSE)$x %>% as.character()
# her2.dn.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.R.output/her2.dn.csv", stringsAsFactors=FALSE)$x %>% as.character()
# basal.up.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.R.output/basal.up.csv", stringsAsFactors=FALSE)$x %>% as.character()
# basal.dn.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.R.output/basal.dn.csv", stringsAsFactors=FALSE)$x %>% as.character()
# 
# M.vs.P.de.gene <- c(lumb.up.gene,lumb.dn.gene,luma.up.gene,luma.dn.gene,her2.up.gene,her2.dn.gene,basal.up.gene,basal.dn.gene) %>% unique
# load('client-side/output/DE.R.output/DE.RData')
# cut.off <- 0.01
# res     <- de.res.liver.vs.breast
# de.res.liver.vs.breast.up.gene <- rownames(res)[res$log2FoldChange > 1  & res$padj < cut.off]
# de.res.liver.vs.breast.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < cut.off]
# 
# pooled.de.gene <- c(de.res.liver.vs.breast.up.gene,de.res.liver.vs.breast.dn.gene) %>% unique
# 
# sum(immune.gene.df$ensemble.id %in% pooled.de.gene)
# 
# liver.vs.breast.immune.signature <- immune.gene.df[immune.gene.df$ensemble.id %in% pooled.de.gene == FALSE,] %>% unique
# liver.vs.breast.immune.signature <- liver.vs.breast.immune.signature[complete.cases(liver.vs.breast.immune.signature),]
