CCLE.breast.cancer.cell.line.characteristic <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/CCLE.breast.cancer.cell.line.characteristic.csv", stringsAsFactors=FALSE)
basal.cell.line                             <- CCLE.breast.cancer.cell.line.characteristic$Cell.line.name[CCLE.breast.cancer.cell.line.characteristic$PAM50.mRNA == 'Basal-like']
#basal.cell.line                             <- paste(basal.cell.line,'BREAST',sep="_")

load('~/Project/Cancer2CellLine/server-side/RData/CCLE_BRACA.RData')

c.cell.line <- intersect(basal.cell.line,colnames(BRACA.log2.fpkm.matrix))
fpkm.matrix <- BRACA.log2.fpkm.matrix[,c.cell.line]
m <- apply(fpkm.matrix,1,median)
fpkm.matrix <- fpkm.matrix[m > 1,]

require(GSVA)
basal.up.gene.immune.excluded <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
basal.dn.gene.immune.excluded <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x

ssgsea.rs <- GSVA::gsva(expr = fpkm.matrix ,gset.idx.list = list(up.gene=basal.up.gene.immune.excluded,dn.gene=basal.dn.gene.immune.excluded), method='ssgsea',ssgsea.norm=FALSE)

up.score <- ssgsea.rs['up.gene',] %>% scale
dn.score <- ssgsea.rs['dn.gene',] %>% scale
