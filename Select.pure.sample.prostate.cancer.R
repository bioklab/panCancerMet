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




cor.cut.off                           <- 0.3
rank.cut.off                          <- 1010


load('server-side/RData/Prostate Adenocarcinoma.RData')
TCGA.prostate.cancer.log2.fpkm.matrix  <- log2.fpkm.matrix
rs.TCGA                                <- pick.out.cell.line(expr.of.samples = TCGA.prostate.cancer.log2.fpkm.matrix,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                               <- rs.TCGA$correlation.matrix[,rs.TCGA$best.cell.line]
cell.line.rank                         <- apply(rs.TCGA$correlation.matrix, 1, function(x) rank(x)[rs.TCGA$best.cell.line] )
pure.PRI.prostate.cancer.sample        <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]



load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
flag                                     <- MET500.sample.meta$cancer.type == 'Prostate Adenocarcinoma' & MET500.sample.meta$LibrarySelection == 'PolyA'  
MET500.sample                            <- MET500.sample.meta$Run[flag]
MET500.prostate.cancer.log2.fpkm.matrix  <- MET500.log2.fpkm.matrix[,MET500.sample]

rs.MET500                                <- pick.out.cell.line(expr.of.samples = MET500.prostate.cancer.log2.fpkm.matrix,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                                 <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                           <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.prostate.cancer.sample          <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]
biopsy.site                              <- MET500.sample.meta[pure.MET.prostate.cancer.sample,'biopsy.site'] %>% as.character()
pure.MET.prostate.cancer.sample          <- pure.MET.prostate.cancer.sample[biopsy.site == 'LIVER']

MET500.best.cell.line <- rs.MET500$best.cell.line
TCGA.best.cell.line   <- rs.TCGA$best.cell.line

save(file = 'client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData',list=c('pure.PRI.prostate.cancer.sample','pure.MET.prostate.cancer.sample','MET500.best.cell.line','TCGA.best.cell.line'))
     
