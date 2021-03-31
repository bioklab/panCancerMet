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



load('server-side/RData/PROSTATE_SRP253428.RData')
liver.met.id         <- PROSTATE_SRP253428_Metadata$Run[PROSTATE_SRP253428_Metadata$source_name == 'CRPC liver metastasis'] %>% as.character()
MET500.prostate.cancer.log2.fpkm.matrix  <- PROSTATE_SRP253428_log2.fpkm.matrix[,liver.met.id]

rs.MET500                                <- pick.out.cell.line(expr.of.samples = MET500.prostate.cancer.log2.fpkm.matrix,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                                 <- rs.MET500$correlation.matrix[,rs.MET500$best.cell.line]
cell.line.rank                           <- apply(rs.MET500$correlation.matrix, 1, function(x) rank(x)[rs.MET500$best.cell.line] )
pure.MET.prostate.cancer.sample.SRP253428  <- names(expr.cor)[(expr.cor > cor.cut.off)  & (cell.line.rank > rank.cut.off) ]


load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')

################################################################################################################
# Prepare primary cancer data
################################################################################################################
load('server-side/RData//Prostate Adenocarcinoma.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.prostate.cancer.sample]
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]


################################################################################################################
# Prepare metastatic cancer data
################################################################################################################
MET.log2.read.count.matrix <- PROSTATE_SRP253428_log2.read.count.matrix[,pure.MET.prostate.cancer.sample.SRP253428]
MET.log2.tpm.matrix        <- PROSTATE_SRP253428_log2.tpm.matrix[,pure.MET.prostate.cancer.sample.SRP253428]


source('client-side/code/DEBoost.R')
load('server-side//RData//Liver.RData')
Male.sample                   <- sample.meta.df$sample.id[sample.meta.df$gender == 'Male'] %>% as.character()
REF.log2.read.count.matrix    <- log2.read.count.matrix[,Male.sample]
REF.log2.tpm.matrix           <- log2.tpm.matrix[,Male.sample]

PRAD.SRP253428.DE.rs <- perform.DE.analysis.between.primary.and.metastatic.cancer(
  PRI.log2.tpm.matrix = PRI.log2.tpm.matrix, PRI.log2.read.count.matrix = PRI.log2.read.count.matrix,
  MET.log2.tpm.matrix = MET.log2.tpm.matrix, MET.log2.read.count.matrix = MET.log2.read.count.matrix,
  REF.log2.tpm.matrix = REF.log2.tpm.matrix, REF.log2.read.count.matrix = REF.log2.read.count.matrix,
  TCGA.best.cell.line = TCGA.best.cell.line, MET500.best.cell.line = rs.MET500$best.cell.line
  
)
save(file='client-side/output/DE.prostate.cancer.SRP253428.R.output/DE.prostate.cancer.SRP253428.RData',list=c('PRAD.SRP253428.DE.rs','pure.MET.prostate.cancer.sample.SRP253428'))
