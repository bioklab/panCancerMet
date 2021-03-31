#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81558  also a good dataset
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110225
require(plyr)
require(dplyr)
require(genefu)
require(Rtsne)
require(ggplot2)
require(dplyr)
source('client-side/code/util.R')


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




load('server-side/RData/COLORECTAL_SRP029880.RData')
met.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'metastatic colorectal cancer to the liver'] %>% as.character()
pri.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'primary colorectal cancer'] %>% as.character()






cor.cut.off <- 0.3
rank.cut.off <- 1010

TCGA.colorectal.cancer.log2.fpkm.matrix  <- COLORECTAL_SRP029880_log2.fpkm.matrix[,pri.sample]
rs.TCGA                                  <- pick.out.cell.line(expr.of.samples = TCGA.colorectal.cancer.log2.fpkm.matrix,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                                 <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                           <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.colorectal.cancer.sample        <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]




MET500.colorectal.cancer.log2.fpkm.matrix  <- COLORECTAL_SRP029880_log2.fpkm.matrix[,met.sample]
rs.MET500                                  <- pick.out.cell.line(expr.of.samples = MET500.colorectal.cancer.log2.fpkm.matrix,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                                   <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                             <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.colorectal.cancer.sample          <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]

MET500.best.cell.line <- rs.MET500$best.cell.line
TCGA.best.cell.line   <- rs.TCGA$best.cell.line


save(file = 'client-side/output/Select.pure.sample.colorectal.cancer.R.output/Select.pure.sample.colorectal.cancer.RData',list=c('pure.PRI.colorectal.cancer.sample','pure.MET.colorectal.cancer.sample','MET500.best.cell.line','TCGA.best.cell.line'))



#######
# met.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'metastatic colorectal cancer to the liver'] %>% as.character()
# pri.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'primary colorectal cancer'] %>% as.character()
# nor.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'normal-looking surrounding colonic epithelium'] %>% as.character()
# 
# 
# pca.rs     <- prcomp(COLORECTAL_SRP029880_log2.fpkm.matrix[,nor.sample] %>% t)
# plot(pca.rs$x[,1:2])
# pc1        <- pca.rs$x[,1]
# pc2         <- pca.rs$x[,2]
# 
# flag <- pc1 < 10 & pc2 > -20
# nor.sample <- (nor.sample)[flag]
# 
# 
# 
# pca.rs        <- prcomp(COLORECTAL_SRP029880_log2.fpkm.matrix[,nor.sample] %>% t)
# plot(pca.rs$x[,1:2])
# 
# pc1.rotation  <- pca.rs$rotation[,'PC1']
# pc2.rotation  <- pca.rs$rotation[,'PC2']
# 
# 
# 
# m.value  <- apply(COLORECTAL_SRP029880_log2.fpkm.matrix[,nor.sample],1,mean)
# 
# met.pc.1 <- apply(COLORECTAL_SRP029880_log2.fpkm.matrix[,pri.sample], 2, function(x)  c((x-m.value) %*% (pc1.rotation)  ))
# met.pc.2 <- apply(COLORECTAL_SRP029880_log2.fpkm.matrix[,pri.sample], 2, function(x)  c((x-m.value) %*% (pc2.rotation)  ))
# points(met.pc.1,met.pc.2,col='purple',pch=19)
# 
# 
# met.pc.1 <- apply(COLORECTAL_SRP029880_log2.fpkm.matrix[,met.sample], 2, function(x)  c((x-m.value) %*% (pc1.rotation)  ))
# met.pc.2 <- apply(COLORECTAL_SRP029880_log2.fpkm.matrix[,met.sample], 2, function(x)  c((x-m.value) %*% (pc2.rotation)  ))
# points(met.pc.1,met.pc.2,col='red',pch=19)
# 
# 
# 
# 
# #pca.rs <- prcomp(COLORECTAL_SRP029880_log2.fpkm.matrix[,c(pri.sample)] %>% t) 
# 
# 
# # x <- (pc1.rotation) %>% sort(decreasing = FALSE) 
# # g <- intersect(names(x)[1:50],median.tpm.matrix %>% rownames)
# # pheatmap(log2(median.tpm.matrix[g,] + 1))
# # 
# # 
# # pc1 <- pca.rs$x[,1]
# # pc2 <- pca.rs$x[,2]
# # 
# 
# # nor.sample <- COLORECTAL_SRP029880_Metadata$Run[COLORECTAL_SRP029880_Metadata$tissue == 'normal-looking surrounding colonic epithelium'] %>% as.character()
# 
# # 
# # 
# # pca.rs <- prcomp(COLORECTAL_SRP029880_log2.fpkm.matrix[,nor.sample] %>% t)
# # plot(pca.rs$x[,1:2])
# # 
# 
# # 
# # 
# # x <- (pc1.rotation) %>% sort(decreasing = TRUE) 
# # g <- intersect(names(x)[1:50],median.tpm.matrix %>% rownames)
# # pheatmap(log2(median.tpm.matrix[g,] + 1))
# 
# 
# # # pc1.rotation.abs <- pca.rs$rotation[,'PC1'] %>% abs %>% sort(decreasing = TRUE)
# # pc1.rotation     <- pca.rs$rotation[,'PC1']
# # pc2.rotation     <- pca.rs$rotation[,'PC2']
# 
# 
# 
# 
# 
# 
# 
# 
# pri.sample <- names(pc1)[pc1 < -20]
# met.sample <- names(met.pc.1)[met.pc.1 < -20]
