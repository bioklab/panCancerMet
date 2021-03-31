require(AnnotationDbi)
library(org.Hs.eg.db)
source('client-side/code/DEBoost.R')

################################################################################################################
# Prepare liver data, here we  want both male and female samples
################################################################################################################
load('server-side//RData//Liver.RData')
REF.log2.tpm.matrix           <- log2.tpm.matrix
REF.log2.read.count.matrix    <- log2.read.count.matrix


load('client-side/output/Select.pure.sample.NET.si.cancer.R.output/Select.pure.sample.NET.si.cancer.RData')
load('server-side/RData/GEP.NET.RData')
################################################################################################################
# Prepare primary cancer data
################################################################################################################
PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.si.cancer.sample]
PRI.log2.tpm.matrix       <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.si.cancer.sample]


################################################################################################################
# Prepare metastatic cancer data
################################################################################################################
MET.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.MET.NET.si.cancer.sample]
MET.log2.tpm.matrix       <- GEP.NET.log2.tpm.matrix[,pure.MET.NET.si.cancer.sample]


NET.SI.DE.rs <- perform.DE.analysis.between.primary.and.metastatic.cancer(
  PRI.log2.tpm.matrix = PRI.log2.tpm.matrix, PRI.log2.read.count.matrix = PRI.log2.read.count.matrix,
  MET.log2.tpm.matrix = MET.log2.tpm.matrix, MET.log2.read.count.matrix = MET.log2.read.count.matrix,
  REF.log2.tpm.matrix = REF.log2.tpm.matrix, REF.log2.read.count.matrix = REF.log2.read.count.matrix,
  TCGA.best.cell.line = TCGA.best.cell.line, MET500.best.cell.line = MET500.best.cell.line
  
)

save(file='client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData',list=c('NET.SI.DE.rs'))


# up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = DE.rs$tumor.intrinsic.DE.gene.rs$up.gene) 
# GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
# GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 & GO.rs.1.up$pvalue < 0.01, ]
# 
# dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene) 
# GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
# GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 & GO.rs.1.dn$pvalue < 0.01, ]
# 



