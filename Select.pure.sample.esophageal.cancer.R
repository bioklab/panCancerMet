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




cor.cut.off <- 0.3

load('server-side/RData/Esophageal Carcinoma.RData')



TCGA.esophageal.cancer.log2.fpkm.matrix  <- log2.fpkm.matrix
pca.rs <- prcomp(TCGA.esophageal.cancer.log2.fpkm.matrix %>% t)
pc1   <- pca.rs$x[,1]
sa    <- names(pc1)[pc1 > 0]

rs.TCGA                                  <- pick.out.cell.line(expr.of.samples = TCGA.esophageal.cancer.log2.fpkm.matrix[,sa],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                                 <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
pure.PRI.esophageal.cancer.sample        <- names(expr.cor)[expr.cor > cor.cut.off]



load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
flag                                       <- MET500.sample.meta$cancer.type == 'Esophageal Adenocarcinoma' & MET500.sample.meta$LibrarySelection == 'PolyA'  
MET500.sample                              <- MET500.sample.meta$Run[flag]
MET500.esophageal.cancer.log2.fpkm.matrix  <- MET500.log2.fpkm.matrix[,MET500.sample]

rs.MET500                                  <- pick.out.cell.line(expr.of.samples = MET500.esophageal.cancer.log2.fpkm.matrix,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                                   <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
pure.MET.esophageal.cancer.sample          <- names(expr.cor)[expr.cor > cor.cut.off]
biopsy.site                                <- MET500.sample.meta[pure.MET.esophageal.cancer.sample,'biopsy.site'] %>% as.character()
pure.MET.esophageal.cancer.sample          <- pure.MET.esophageal.cancer.sample[biopsy.site == 'LIVER']


save(file = 'client-side/output/Select.pure.sample.esophageal.cancer.output/Select.pure.sample.esophageal.cancer.RData',list=c('pure.PRI.esophageal.cancer.sample','pure.MET.esophageal.cancer.sample'))

