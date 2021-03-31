##### Aim: compute TCGA.sample - CCLE.cell.line correlation values, as an estimator of tumor purity

require(plyr)
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


load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix
TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix


MET500.sample <- MET500.breast.cancer.polyA.Basal.sample
TCGA.sample   <- TCGA.breast.cancer.polyA.Basal.sample
rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
purity.MET500 <- rs.MET500$correlation.matrix[,'HCC70_BREAST'] # use HCC70 cell line
purity.TCGA   <- rs.TCGA$correlation.matrix[,'HCC1569_BREAST']
#purity.MET500 <- rs.MET500$correlation.matrix[,'HCC1569_BREAST'] # use HCC70 cell line
#purity.TCGA   <- rs.TCGA$correlation.matrix[,'HCC1569_BREAST']
tumor.purity.based.on.cell.line.vec <- c(purity.MET500,purity.TCGA)


MET500.sample <- MET500.breast.cancer.polyA.Her2.sample
TCGA.sample   <- TCGA.breast.cancer.polyA.Her2.sample
rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
purity.MET500 <- rs.MET500$correlation.matrix[,'EFM192A_BREAST'] # use EFM192A cell line
purity.TCGA   <- rs.TCGA$correlation.matrix[,'EFM192A_BREAST']
tumor.purity.based.on.cell.line.vec <- c(tumor.purity.based.on.cell.line.vec,purity.MET500,purity.TCGA)



MET500.sample <- MET500.breast.cancer.polyA.LumA.sample
TCGA.sample   <- TCGA.breast.cancer.polyA.LumA.sample
rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
purity.MET500 <- rs.MET500$correlation.matrix[,'MDAMB415_BREAST'] # use MDAMB415 cell line
purity.TCGA   <- rs.TCGA$correlation.matrix[,'MDAMB415_BREAST']
tumor.purity.based.on.cell.line.vec <- c(tumor.purity.based.on.cell.line.vec,purity.MET500,purity.TCGA)



MET500.sample <- MET500.breast.cancer.polyA.LumB.sample
TCGA.sample   <- TCGA.breast.cancer.polyA.LumB.sample
rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
purity.MET500 <- rs.MET500$correlation.matrix[,'BT483_BREAST'] # use BT483 cell line
purity.TCGA   <- rs.TCGA$correlation.matrix[,'BT483_BREAST']
tumor.purity.based.on.cell.line.vec <- c(tumor.purity.based.on.cell.line.vec,purity.MET500,purity.TCGA)

save(list = 'tumor.purity.based.on.cell.line.vec',file='client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')




