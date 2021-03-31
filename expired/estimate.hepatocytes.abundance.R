require(plyr)
require(dplyr)
require(foreach)


#KDM5D <- 'ENSG00000012817'  Y chrosome linked gene
#KDM5D.expr                  <- GTex.liver.log2.fpkm.matrix[KDM5D,]

###################################
load('server-side//RData//Liver.RData')
GTex.liver.log2.fpkm.matrix <- log2.fpkm.matrix
female.sample               <- sample.meta.df$sample.id[sample.meta.df$gender == 'Female'] %>% as.character()
GTex.liver.log2.fpkm.matrix <- GTex.liver.log2.fpkm.matrix[,female.sample]
m.expr                      <- apply(GTex.liver.log2.fpkm.matrix,1,median) 
liver.expressed.gene        <- names(m.expr)[m.expr > log2(10+1)]

load('server-side/RData//Breast Invasive Carcinoma.RData')
TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
m.expr                                    <- apply(TCGA.breast.cancer.log2.fpkm.matrix,1,median) 
bc.not.expressed.gene                     <- names(m.expr)[m.expr < log2(0.1+1)]

liver.specific.gene <- intersect(liver.expressed.gene,bc.not.expressed.gene)

liver.expr <- apply(GTex.liver.log2.fpkm.matrix[liver.specific.gene,],1,median) 


load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 

ME500.brca.sample <- intersect(MET500.liver.sample,c(MET500.breast.cancer.polyA.Basal.sample,MET500.breast.cancer.polyA.Her2.sample,MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.LumB.sample))

hepatocyte.abundance.vec <- foreach(s = ME500.brca.sample,.combine='c') %do% {
    met.expr        <- MET500.log2.fpkm.matrix[liver.specific.gene,s]
    df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
    flag            <- df$liver.expr >0 & df$met.expr > 0
    df              <- df[flag,]
    df$met.expr     <- log2(2^df$met.expr - 1)
    df$liver.expr   <- log2(2^df$liver.expr - 1)
    intercept       <- median(df$met.expr - df$liver.expr )
    2^intercept
}
names(hepatocyte.abundance.vec) <- ME500.brca.sample

hepatocyte.abundance.df <- foreach(s = ME500.brca.sample) %do% {
    met.expr        <- MET500.log2.fpkm.matrix[liver.specific.gene,s]
    df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
    flag            <- df$liver.expr >0 & df$met.expr > 0
    df              <- df[flag,]
    df$met.expr     <- log2(2^df$met.expr - 1)
    df$liver.expr   <- log2(2^df$liver.expr - 1)
    df
}
names(hepatocyte.abundance.df) <- ME500.brca.sample
save(file='client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData',list = c('hepatocyte.abundance.vec','hepatocyte.abundance.df'))
