require(plyr)
require(dplyr)
require(foreach)
require(AnnotationDbi)
library(org.Hs.eg.db)

############################ select liver-expressed gene ####################
load('server-side/RData/SRP068976_PairNor.RData')
load('server-side/RData/SRP174991_PairNor.RData')

liver.expr.cut.off <- 2

liver.log2.fpkm.matrix           <- SRP174991_log2.fpkm.matrix
KDM5D                            <- 'ENSG00000012817'
KDM5D.expr                       <- liver.log2.fpkm.matrix[KDM5D,]
plot(KDM5D.expr %>% sort)
SRP174991.male.sample            <- names(KDM5D.expr)[KDM5D.expr > log2(2+1)]
SRP174991.female.sample          <- names(KDM5D.expr)[KDM5D.expr < log2(2+1)]
female.m.expr                    <- apply(liver.log2.fpkm.matrix[,SRP174991.female.sample],1,median) 
SRP174991.liver.expressed.gene   <- names(female.m.expr)[female.m.expr > log2(liver.expr.cut.off+1)]

liver.log2.fpkm.matrix           <- SRP068976_log2.fpkm.matrix
KDM5D                            <- 'ENSG00000012817'
KDM5D.expr                       <- liver.log2.fpkm.matrix[KDM5D,]
plot(KDM5D.expr %>% sort)
SRP068976.male.sample            <- names(KDM5D.expr)[KDM5D.expr > log2(2+1)]
SRP068976.female.sample          <- names(KDM5D.expr)[KDM5D.expr < log2(2+1)]
female.m.expr                    <- apply(liver.log2.fpkm.matrix[,SRP068976.female.sample],1,median) 
SRP068976.liver.expressed.gene   <- names(female.m.expr)[female.m.expr > log2(liver.expr.cut.off+1)]


load('server-side//RData//Liver.RData')
GTex.log2.fpkm.matrix            <- log2.fpkm.matrix
liver.log2.fpkm.matrix           <- GTex.log2.fpkm.matrix
KDM5D                            <- 'ENSG00000012817'
KDM5D.expr                       <- liver.log2.fpkm.matrix[KDM5D,]
plot(KDM5D.expr %>% sort)
GTex.male.sample            <- names(KDM5D.expr)[KDM5D.expr > log2(2+1)]
GTex.female.sample          <- names(KDM5D.expr)[KDM5D.expr < log2(2+1)]
female.m.expr               <- apply(liver.log2.fpkm.matrix[,GTex.female.sample],1,median) 
GTex.liver.expressed.gene   <- names(female.m.expr)[female.m.expr > log2(liver.expr.cut.off+1)]



liver.expressed.gene <- intersect(SRP174991.liver.expressed.gene,SRP068976.liver.expressed.gene)
#liver.expressed.gene <- GTex.liver.expressed.gene

######################select genes NOT expressed in brca ###################
load('server-side/RData//Breast Invasive Carcinoma.RData')
TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
m.expr                                    <- apply(TCGA.breast.cancer.log2.fpkm.matrix,1,median) 
bc.not.expressed.gene                     <- names(m.expr)[m.expr < log2(0.01+1)]


################# get the liver-specific expressed gene ###############
liver.specific.gene    <- intersect(liver.expressed.gene,bc.not.expressed.gene)
protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
liver.specific.gene    <- intersect(liver.specific.gene,protein.coding.gene.id)

#liver.expr             <- apply(SRP174991_log2.fpkm.matrix[liver.specific.gene,SRP174991.female.sample],1,median)
liver.expr             <- apply(SRP068976_log2.fpkm.matrix[liver.specific.gene,SRP068976.female.sample],1,median)
#liver.expr             <- apply(GTex.log2.fpkm.matrix[liver.specific.gene,GTex.female.sample],1,median)


##################### organize breast cancer liver metastsis data from GEO #############
load('server-side/RData/BRACA_SRP043470.RData')
SRP043470.met.sample      <- BRACA_SRP043470_Metadata$Run[BRACA_SRP043470_Metadata$source_name=='Liver Metastasis Tumour'] %>% as.character()

load('server-side/RData/BRACA_LIVERMET.RData')
LIVERMET.log2.fpkm.matrix <- log2.fpkm.matrix
SRP102119.met.sample      <- c('SRR5357752','SRR5357760','SRR5357767') #From project SRP102119
SRP125109.met.sample      <- c('SRR6298032','SRR6298035','SRR6298037','SRR6298040','SRR6298041') # From project SRP125109

METGEO.log2.fpkm.matrix <- cbind(BRACA_SRP043470_log2.fpkm.matrix,LIVERMET.log2.fpkm.matrix[BRACA_SRP043470_log2.fpkm.matrix %>% rownames,])
METGEO.log2.fpkm.matrix <- METGEO.log2.fpkm.matrix[,c(SRP043470.met.sample,SRP102119.met.sample,SRP125109.met.sample)]
GEO.brca.sample         <- colnames(METGEO.log2.fpkm.matrix)

####################### let us perform the analysis ####################
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
ME500.brca.sample    <- intersect(MET500.liver.sample,c(MET500.breast.cancer.polyA.Basal.sample,MET500.breast.cancer.polyA.Her2.sample,MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.LumB.sample))

MET.log2.fpkm.matrix <- cbind(METGEO.log2.fpkm.matrix,MET500.log2.fpkm.matrix[rownames(METGEO.log2.fpkm.matrix),ME500.brca.sample])


MET.rs <- foreach(s = c(ME500.brca.sample,GEO.brca.sample),.combine='rbind') %do% {
    df              <- data.frame(liver.expr = liver.expr[liver.specific.gene],met.expr=MET.log2.fpkm.matrix[liver.specific.gene,s])
    flag            <- df$liver.expr >0 & df$met.expr > 0
    df              <- df[flag,]
    df$met.expr     <- log2(2^df$met.expr - 1)
    df$liver.expr   <- log2(2^df$liver.expr - 1)
    intercept       <- median(df$met.expr - df$liver.expr )
    df$residual     <- df$met.expr - (df$liver.expr + intercept)
    df$gene         <- rownames(df)
    df$sample       <- s
    sd.estimate     <- mad(df$met.expr - df$liver.expr) *  1.4826
    df$dn.p.value   <- pnorm(df$residual,mean=0,sd=sd.estimate,lower.tail = TRUE)
    df$up.p.value   <- pnorm(df$residual,mean=0,sd=sd.estimate,lower.tail = FALSE)
    df$intercept    <- intercept
    df$sd.estimate  <- sd.estimate
    df$cor          <- cor(df$met.expr,df$liver.expr,method='spearman')
    df[,c('sample','intercept','sd.estimate','cor','gene','met.expr','liver.expr','residual','up.p.value','dn.p.value')]
}

combine.p.value <- function(x) {
    # q     <- quantile(x$residual)
    # lower <- q[2] - 1.5 * IQR(x$residual)
    # upper <- q[4] + 1.5 * IQR(x$residual)
    # x     <- x[x$residual > lower & x$residual < upper,]
    chi.sequare.statistic <- -2 * sum(log(x$up.p.value) )
    up.p.value <- pchisq(chi.sequare.statistic,df = 2 * nrow(x),lower.tail=FALSE)
    chi.sequare.statistic <- -2 * sum(log(x$dn.p.value) )
    dn.p.value <- pchisq(chi.sequare.statistic,df = 2 * nrow(x),lower.tail=FALSE)
    data.frame(up.p.value=up.p.value,dn.p.value=dn.p.value)
}


MET.rs                   <- MET.rs[MET.rs$cor > 0.8 & MET.rs$intercept < -1 ,]
MET.p.value.df           <- ddply(MET.rs,.(gene),combine.p.value)
rownames(MET.p.value.df) <- MET.p.value.df$gene
MET.p.value.df$up.fdr    <- p.adjust(MET.p.value.df$up.p.value,method='fdr')
MET.p.value.df$dn.fdr    <- p.adjust(MET.p.value.df$dn.p.value,method='fdr')


MET.up.gene    <- MET.p.value.df[MET.p.value.df$up.fdr < 0.05,'gene']
MET.dn.gene    <- MET.p.value.df[MET.p.value.df$dn.fdr < 0.05,'gene']
MET.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = MET.up.gene)
MET.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = MET.dn.gene)

# GEO.up.gene     <- GEO.p.value.df[GEO.p.value.df$up.p.value   < 0.05,'gene']
# GEO.dn.gene     <- GEO.p.value.df[GEO.p.value.df$dn.p.value   < 0.05,'gene']
# common.up.gene  <- intersect(MET500.up.gene,GEO.up.gene)
# common.dn.gene  <- intersect(MET500.dn.gene,GEO.dn.gene)


# MET500.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = MET500.up.gene)
# MET500.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = MET500.dn.gene)
# GEO.up.gene.df    <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = GEO.up.gene)
# GEO.dn.gene.df    <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = GEO.dn.gene)
# common.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = common.up.gene)
# common.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = common.dn.gene)

###########
DE between samples with high and low cor.value

#################
flag       <- MET500.p.value.df.SRP068976$up.fdr < 0.05 &  MET500.p.value.df.SRP174991$up.fdr < 0.05
up.gene    <- MET500.p.value.df.SRP068976$gene[flag]
MET500.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)

flag       <- MET500.p.value.df.SRP068976$dn.fdr < 0.05 &  MET500.p.value.df.SRP174991$dn.fdr < 0.05
dn.gene    <- MET500.p.value.df.SRP068976$gene[flag]
MET500.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = dn.gene)



flag       <- GEO.p.value.df.SRP068976$up.p.value < 0.05 &  GEO.p.value.df.SRP174991$up.p.value < 0.05
up.gene    <- GEO.p.value.df.SRP068976$gene[flag]
GEO.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)

flag       <- GEO.p.value.df.SRP068976$dn.p.value < 0.05 &  GEO.p.value.df.SRP174991$dn.p.value < 0.05
dn.gene    <- GEO.p.value.df.SRP068976$gene[flag]
GEO.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = dn.gene)






require(idr)
x <- cbind(-1 * (MET500.p.value.df$up.p.value %>% log10), -1 * (GEO.p.value.df$up.p.value %>% log10))
uv <- get.correspondence(x[,1],x[,2],seq(0.01, 0.99, by=1/28))

plot(uv$psi$t, uv$psi$value, xlab="t", ylab="psi", xlim=c(0, max(uv$psi$t)),
     ylim=c(0, max(uv$psi$value)), cex.lab=2)
lines(uv$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)

# rownames(x) <- MET500.p.value.df$gene
# mu <- 2.6
# sigma <- 1.3
# rho <- 0.8
# p <- 0.7
# idr.out <- est.IDR(x, mu, sigma, rho, p, eps=0.001, max.ite=20)
# 
# IDR.level <- 0.00001
# x.selected <- select.IDR(x, idr.out$IDR, IDR.level)
#################################### Trash code below ######################################
# g                 <- intersect(MET500.p.value.df$gene,GEO.p.value.df$gene)
# MET500.p.value.df <- MET500.p.value.df[g,]
# GEO.p.value.df    <- GEO.p.value.df[g,]
# plot(x=MET500.p.value.df$up.p.value %>% log10() %>% rank,y=GEO.p.value.df$up.p.value %>% log10() %>% rank)
# plot(x=MET500.p.value.df$dn.p.value %>% log10() %>% rank,y=GEO.p.value.df$dn.p.value %>% log10() %>% rank)
# 
# # 
# 
# KDM5D <- 'ENSG00000012817'
# CYP1A1 <- 'ENSG00000140465'
# KDM5D.expr <- SRP174991_log2.fpkm.matrix[KDM5D,]
# male.sample <- names(KDM5D.expr)[KDM5D.expr >2]
# female.sample <- names(KDM5D.expr)[KDM5D.expr < 2]
# 
# tmp1 <- apply(SRP174991_log2.fpkm.matrix[liver.specific.gene,female.sample],1,median)
# tmp2 <- apply(SRP174991_log2.fpkm.matrix[liver.specific.gene,male.sample],1,median)
# 
# liver.specific.gene <- names(tmp1)[tmp1 > log2(10 +1) ]
# #################################
# require(dplyr)
# ##############################  select liver-specific genes ########################
require(CePa)
tissue.median.tpm.matrix <- read.gct(file = 'client-side/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')

get.gene.id <- function(x){
  strsplit(x,split='\\.')[[1]] %>% unlist  %>% head(1)
}
rownames(tissue.median.tpm.matrix) <- sapply(rownames(tissue.median.tpm.matrix),get.gene.id)






is.liver.specific <- function(x) {
    liver.expr <- x['Liver']
    y          <- sort(x,decreasing = TRUE)
    if((liver.expr < 1) | (liver.expr != max(y))){
        return(FALSE)
    }else{
        if(y[2] ==0 | liver.expr/y[2] > 5){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }

}
flag                       <- apply(tissue.median.tpm.matrix,1,is.liver.specific)
liver.specific.gene.across.tissues        <- rownames(tissue.median.tpm.matrix)[flag]
names(liver.specific.gene.across.tissues) <- NULL



# 
# load('server-side//RData//Liver.RData')
# GTex.liver.log2.fpkm.matrix       <- log2.fpkm.matrix
# liver.specific.gene               <- intersect(liver.specific.gene,rownames(GTex.liver.log2.fpkm.matrix))
# liver.expr                        <- apply(GTex.liver.log2.fpkm.matrix[liver.specific.gene,],1,median)
# liver.specific.gene               <- names(liver.expr)[liver.expr > log2(5+1)]
# 



# SAA1 <- 'ENSG00000173432'
# DHODH <- 'ENSG00000102967'
# 
# g <- 'ENSG00000276076'
# MET500.rs[MET500.rs$sample %in% MET500.fit.sample & MET500.rs$gene == g,]
# mad.df <- ddply(MET500.rs,.(gene),function(x) mad(x$residual))
# mad.df <- mad.df[order(mad.df$V1,decreasing = TRUE),]
# # 
# # median.df <- ddply(rs,.(gene),function(x) median(x$residual))
# # median.df <- median.df[order(median.df$V1,decreasing = TRUE),]
# 
# 
# SAA1 <- 'ENSG00000173432'
# DHODH <- 'ENSG00000102967'
# 
# draw.df <- GEO.rs[GEO.rs$gene == SAA1,]
# draw.df$subtype <- 'HE'
# #
# draw.df[draw.df$sample %in% MET500.breast.cancer.polyA.Basal.sample,'subtype'] <- 'Basal'
# draw.df[draw.df$sample %in% MET500.breast.cancer.polyA.LumA.sample,'subtype']  <- 'LumA'
# draw.df[draw.df$sample %in% MET500.breast.cancer.polyA.LumB.sample,'subtype']  <- 'LumB'
# draw.df[draw.df$sample %in% MET500.breast.cancer.polyA.Her2.sample,'subtype']  <- 'Her2'
# #
# draw.df <- draw.df[order(draw.df$residual),]
# draw.df$rank <- 1:nrow(draw.df)
# ggplot(draw.df[order(draw.df$residual),]) + geom_point(aes(y=residual,x=rank,color=subtype))
# 
# 
# 
# 
# # ################################################################3
# # MET500.sample   <- MET500.breast.cancer.polyA.LumA.sample
# # MET500.sample   <- intersect(MET500.sample,MET500.liver.sample)
# # met.expr        <- apply(MET500.log2.fpkm.matrix[liver.specific.gene,MET500.sample],1,median)
# # df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
# # flag            <- df$liver.expr >0 & df$met.expr > 0
# # df              <- df[flag,]
# # df$met.expr     <- log2(2^df$met.expr - 1)
# # df$liver.expr   <- log2(2^df$liver.expr - 1)
# # intercept       <-  median(df$met.expr - df$liver.expr )
# # df$residual     <- df$met.expr - (df$liver.expr + intercept)
# # ggplot(df,aes(x=liver.expr,y=met.expr)) + geom_point() + geom_abline(slope=1,intercept = intercept) 
# # df.luma         <- df
# # 
# # 
# # MET500.sample   <- MET500.breast.cancer.polyA.LumB.sample
# # MET500.sample   <- intersect(MET500.sample,MET500.liver.sample)
# # met.expr        <- apply(MET500.log2.fpkm.matrix[liver.specific.gene,MET500.sample],1,median)
# # df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
# # flag            <- df$liver.expr >0 & df$met.expr > 0
# # df              <- df[flag,]
# # df$met.expr     <- log2(2^df$met.expr - 1)
# # df$liver.expr   <- log2(2^df$liver.expr - 1)
# # intercept       <-  median(df$met.expr - df$liver.expr )
# # df$residual     <- df$met.expr - (df$liver.expr + intercept)
# # ggplot(df,aes(x=liver.expr,y=met.expr)) + geom_point() + geom_abline(slope=1,intercept = intercept) 
# # df.lumb         <- df
# # 
# # 
# # 
# # MET500.sample   <- MET500.breast.cancer.polyA.Basal.sample
# # MET500.sample   <- intersect(MET500.sample,MET500.liver.sample)
# # met.expr        <- apply(MET500.log2.fpkm.matrix[liver.specific.gene,MET500.sample],1,median)
# # df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
# # flag            <- df$liver.expr >0 & df$met.expr > 0
# # df              <- df[flag,]
# # df$met.expr     <- log2(2^df$met.expr - 1)
# # df$liver.expr   <- log2(2^df$liver.expr - 1)
# # intercept       <-  median(df$met.expr - df$liver.expr )
# # df$residual     <- df$met.expr - (df$liver.expr + intercept)
# # ggplot(df,aes(x=liver.expr,y=met.expr)) + geom_point() + geom_abline(slope=1,intercept = intercept) 
# # df.basal         <- df
# # 
# # 
# # MET500.sample   <- MET500.breast.cancer.polyA.Her2.sample
# # MET500.sample   <- intersect(MET500.sample,MET500.liver.sample)
# # met.expr        <- apply(MET500.log2.fpkm.matrix[liver.specific.gene,MET500.sample],1,median)
# # df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
# # flag            <- df$liver.expr >0 & df$met.expr > 0
# # df              <- df[flag,]
# # df$met.expr     <- log2(2^df$met.expr - 1)
# # df$liver.expr   <- log2(2^df$liver.expr - 1)
# # intercept       <-  median(df$met.expr - df$liver.expr )
# # df$residual     <- df$met.expr - (df$liver.expr + intercept)
# # ggplot(df,aes(x=liver.expr,y=met.expr)) + geom_point() + geom_abline(slope=1,intercept = intercept) 
# # df.her2         <- df
# # 
# # 
# # 
# # 
# # 
# # MET500.sample   <- c(MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.LumB.sample)
# # MET500.sample   <- intersect(MET500.sample,MET500.liver.sample)
# # met.expr        <- apply(MET500.log2.fpkm.matrix[liver.specific.gene,MET500.sample],1,median)
# # df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
# # flag            <- df$liver.expr >0 & df$met.expr > 0
# # df              <- df[flag,]
# # df$met.expr     <- log2(2^df$met.expr - 1)
# # df$liver.expr   <- log2(2^df$liver.expr - 1)
# # intercept       <-  median(df$met.expr - df$liver.expr )
# # df$residual     <- df$met.expr - (df$liver.expr + intercept)
# # ggplot(df,aes(x=liver.expr,y=met.expr)) + geom_point() + geom_abline(slope=1,intercept = intercept) 
# # df.pool        <- df
# 
# 
# 
# 
# 
# 
# 
# load('server-side/RData/BRACA_SRP043470.RData')
# m.sample         <- BRACA_SRP043470_Metadata$Run[BRACA_SRP043470_Metadata$source_name=='Liver Metastasis Tumour'] %>% as.character()
# #SRP043470.m.expr <- apply(BRACA_SRP043470_log2.fpkm.matrix[liver.specific.gene,m.sample],1,median)
# rs <- foreach(s = m.sample,.combine='rbind') %do% {
#     met.expr        <- BRACA_SRP043470_log2.fpkm.matrix[liver.specific.gene,s]
#     df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
#     flag            <- df$liver.expr >0 & df$met.expr > 0
#     df              <- df[flag,]
#     df$met.expr     <- log2(2^df$met.expr - 1)
#     df$liver.expr   <- log2(2^df$liver.expr - 1)
#     intercept       <- median(df$met.expr - df$liver.expr )
#     df$residual     <- df$met.expr - (df$liver.expr + intercept)
#     ggplot(df,aes(x=liver.expr,y=met.expr)) + geom_point() + geom_abline(slope=1,intercept = intercept) 
#     df$gene         <- rownames(df)
#     df$sample       <- s
#     sd.estimate     <- mad(df$residual) *  1.4826
#     df$dn.p.value   <- pnorm(df$residual,mean=0,sd=sd.estimate,lower.tail = TRUE)
#     df$up.p.value   <- pnorm(df$residual,mean=0,sd=sd.estimate,lower.tail = FALSE)
#     df
# }
# SRP043470.p.value.df <- ddply(rs,.(gene),combine.p.value)
# rownames(SRP043470.p.value.df) <- SRP043470.p.value.df$gene
# 
# 
# 
# 
# load('server-side/RData/BRACA_LIVERMET.RData')
# LIVERMET.log2.fpkm.matrix <- log2.fpkm.matrix
# m.sample <- c('SRR5357752','SRR5357760','SRR5357767') #From project SRP102119
# #m.sample <- c('SRR6298032','SRR6298035','SRR6298037','SRR6298040','SRR6298041')
# 
# rs <- foreach(s = m.sample,.combine='rbind') %do% {
#   met.expr        <- LIVERMET.log2.fpkm.matrix[liver.specific.gene,s]
#   df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
#   flag            <- df$liver.expr >0 & df$met.expr > 0
#   df              <- df[flag,]
#   df$met.expr     <- log2(2^df$met.expr - 1)
#   df$liver.expr   <- log2(2^df$liver.expr - 1)
#   intercept       <- median(df$met.expr - df$liver.expr )
#   df$residual     <- df$met.expr - (df$liver.expr + intercept)
#   ggplot(df,aes(x=liver.expr,y=met.expr)) + geom_point() + geom_abline(slope=1,intercept = intercept) 
#   df$gene         <- rownames(df)
#   df$sample       <- s
#   sd.estimate     <- mad(df$residual) *  1.4826
#   df$dn.p.value   <- pnorm(df$residual,mean=0,sd=sd.estimate,lower.tail = TRUE)
#   df$up.p.value   <- pnorm(df$residual,mean=0,sd=sd.estimate,lower.tail = FALSE)
#   df
# }
# SRP102119.p.value.df           <- ddply(rs,.(gene),combine.p.value)
# rownames(SRP102119.p.value.df) <- SRP102119.p.value.df$gene
# 
# 
# 
# m.sample <- c('SRR6298032','SRR6298035','SRR6298037','SRR6298040','SRR6298041') # TNBC
# rs <- foreach(s = m.sample,.combine='rbind') %do% {
#   met.expr        <- LIVERMET.log2.fpkm.matrix[liver.specific.gene,s]
#   df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
#   flag            <- df$liver.expr >0 & df$met.expr > 0
#   df              <- df[flag,]
#   df$met.expr     <- log2(2^df$met.expr - 1)
#   df$liver.expr   <- log2(2^df$liver.expr - 1)
#   intercept       <- median(df$met.expr - df$liver.expr )
#   df$residual     <- df$met.expr - (df$liver.expr + intercept)
#   ggplot(df,aes(x=liver.expr,y=met.expr)) + geom_point() + geom_abline(slope=1,intercept = intercept) 
#   df$gene         <- rownames(df)
#   df$sample       <- s
#   sd.estimate     <- mad(df$residual) *  1.4826
#   df$dn.p.value   <- pnorm(df$residual,mean=0,sd=sd.estimate,lower.tail = TRUE)
#   df$up.p.value   <- pnorm(df$residual,mean=0,sd=sd.estimate,lower.tail = FALSE)
#   df
# }
# SRP125109.p.value.df           <- ddply(rs,.(gene),combine.p.value)
# rownames(SRP125109.p.value.df)<-  SRP125109.p.value.df$gene
# 
# 
# 
# 
# 
# 
# 
# get.outlier.gene <- function(x) {
#     iqr <- IQR(x$residual)  
#     q1  <- quantile(x$residual)[2]
#     q3  <- quantile(x$residual)[4]
#     up  <- q3 + 1.5*iqr
#     dn  <- q1 - 1.5*iqr
#     list(up.gene=rownames(x)[x$residual > up],dn.gene=rownames(x)[x$residual < dn])
# }
# 
# luma.outlier.rs      <- get.outlier.gene(df.luma)
# lumb.outlier.rs      <-  get.outlier.gene(df.lumb)
# her2.outlier.rs      <-  get.outlier.gene(df.her2)
# basal.outlier.rs     <-  get.outlier.gene(df.basal)
# SRP043470.outlier.rs <-  get.outlier.gene(df.SRP043470)
# SRP102119.outlier.rs <-  get.outlier.gene(df.SRP102119)
# SRP125109.outlier.rs <-  get.outlier.gene(df.SRP125109)
# 
# up.gene.freq <- c(luma.outlier.rs$up.gene,lumb.outlier.rs$up.gene,her2.outlier.rs$up.gene,basal.outlier.rs$up.gene,SRP043470.outlier.rs$up.gene,SRP102119.outlier.rs$up.gene) %>% table %>% as.data.frame
# dn.gene.freq <- c(luma.outlier.rs$dn.gene,lumb.outlier.rs$dn.gene,her2.outlier.rs$dn.gene,basal.outlier.rs$dn.gene,SRP043470.outlier.rs$dn.gene,SRP102119.outlier.rs$dn.gene) %>% table %>% as.data.frame
# 
# 
# 
# #############
# load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
# 
# 
# g  <- intersect(rownames(de.res.liver.vs.breast.luma),rownames(de.res.metastasis.liver.vs.breast.luma))
# df <- data.frame(x=de.res.liver.vs.breast.luma[g,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.luma[g,'log2FoldChange'])
# rownames(df) <- g
# df <- df[df$x >=5,]
# ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5)  + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_smooth(col='red',size=2)
# 
# loess.fit         <- loess(data = df,formula=y~x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# adjusted.fc       <- loess.fit$residuals
# z.score           <- (adjusted.fc - median(adjusted.fc))/(1.4826 * mad(adjusted.fc))
# 
# luma.dn.gene <- names(adjusted.fc)[adjusted.fc <= -1 & z.score <= -3]
# dn.gene.annotation.df <-  AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = luma.dn.gene)
# 
# 
# 
# 
# g  <- intersect(rownames(de.res.liver.vs.breast.lumb),rownames(de.res.metastasis.liver.vs.breast.lumb))
# df <- data.frame(x=de.res.liver.vs.breast.lumb[g,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.lumb[g,'log2FoldChange'])
# rownames(df) <- g
# df <- df[df$x >=5,]
# ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5)  + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_smooth(col='red',size=2)
# 
# loess.fit         <- loess(data = df,formula=y~x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# adjusted.fc       <- loess.fit$residuals
# z.score           <- (adjusted.fc - median(adjusted.fc))/(1.4826 * mad(adjusted.fc))
# 
# lumb.up.gene <- names(adjusted.fc)[adjusted.fc >=1 & z.score >= 3]
# up.gene.annotation.df <-  AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)
# 
# lumb.dn.gene <- names(adjusted.fc)[adjusted.fc <= -1 & z.score <= -3]
# dn.gene.annotation.df <-  AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = lumb.dn.gene)
# 
# 
# 
# 
# g  <- intersect(rownames(de.res.liver.vs.breast.her2),rownames(de.res.metastasis.liver.vs.breast.her2))
# df <- data.frame(x=de.res.liver.vs.breast.her2[g,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.her2[g,'log2FoldChange'])
# rownames(df) <- g
# df <- df[df$x >=5,]
# ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5)  + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_smooth(col='red',size=2)
# 
# loess.fit         <- loess(data = df,formula=y~x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# adjusted.fc       <- loess.fit$residuals
# z.score           <- (adjusted.fc - median(adjusted.fc))/(1.4826 * mad(adjusted.fc))
# 
# her2.up.gene <- names(adjusted.fc)[adjusted.fc >=1 & z.score >= 3]
# up.gene.annotation.df <-  AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)
# her2.dn.gene <- names(adjusted.fc)[adjusted.fc <= -1 & z.score <= -3]
# dn.gene.annotation.df <-  AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = her2.dn.gene)
# 
# 
# 
# 
# 
# 
# g  <- intersect(rownames(de.res.liver.vs.breast.basal),rownames(de.res.metastasis.liver.vs.breast.basal))
# df <- data.frame(x=de.res.liver.vs.breast.basal[g,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.basal[g,'log2FoldChange'])
# rownames(df) <- g
# df <- df[df$x >=5,]
# ggplot(df,aes(x=x,y=y)) + geom_point(size=2.5)  + ylab('log2FC (MET.vs.PRI)') + xlab('log2FC (LIVER.vs.PRI)') + xlim(c(-15,22)) + ylim(c(-15,22)) + geom_abline(intercept = 0,slope=1) + geom_smooth(col='red',size=2)
# 
# loess.fit         <- loess(data = df,formula=y~x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# adjusted.fc       <- loess.fit$residuals
# z.score           <- (adjusted.fc - median(adjusted.fc))/(1.4826 * mad(adjusted.fc))
# 
# basal.up.gene <- names(adjusted.fc)[adjusted.fc >=1 & z.score >= 3]
# up.gene.annotation.df <-  AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)
# basal.dn.gene <- names(adjusted.fc)[adjusted.fc <= -1 & z.score <= -3]
# dn.gene.annotation.df <-  AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = basal.dn.gene)
# 
# 
# 
# 
# pool.gene <- c(luma.up.gene,lumb.up.gene,her2.up.gene,basal.up.gene)
# freq.df <- table(pool.gene) %>% as.data.frame
# 
# ####################
# require(DESeq2)
# require(dplyr)
# library (VennDiagram)
# require(gplots)
# require(foreach)
# source('client-side/code/util.R')
# require(AnnotationDbi)
# library(org.Hs.eg.db)
# protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
# 
# 
# 
# load('server-side//RData//Liver.RData')
# GTex.liver.log2.read.count.matrix <- log2.read.count.matrix
# GTex.liver.log2.fpkm.matrix       <- log2.fpkm.matrix
# liver.expressed.gene              <- get.expressed.gene(GTex.liver.log2.fpkm.matrix)
# 
# # load('server-side//RData//Breast - Mammary Tissue.RData')
# # GTex.breast.log2.read.count.matrix <- log2.read.count.matrix
# # GTex.breast.log2.fpkm.matrix       <- log2.fpkm.matrix
# # breast.expressed.gene              <- get.expressed.gene(GTex.breast.log2.fpkm.matrix)
# 
# 
# 
# load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
# load('server-side/RData//Breast Invasive Carcinoma.RData')
# load('client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')
# TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix
# TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
# TCGA.sample                               <- pure.TCGA.breast.cancer.polyA.LumB.sample
# TCGA.expressed.gene                       <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample])
# 
# 
# load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
# load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
# MET500.liver.sample     <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# MET500.sample           <- MET500.breast.cancer.polyA.LumB.sample
# MET500.sample           <- intersect(MET500.sample,MET500.liver.sample)
# MET500.expressed.gene   <- get.expressed.gene(MET500.log2.fpkm.matrix[,MET500.sample])
# 
# 
# 
# 
# 
# 
# 
# 
# expressed.gene <- intersect(c(liver.expressed.gene,MET500.expressed.gene) %>% unique,protein.coding.gene.id)
# expr.matrix    <- cbind(MET500.log2.read.count.matrix[expressed.gene,MET500.sample],GTex.liver.log2.read.count.matrix[expressed.gene,])
# expr.matrix    <- 2^expr.matrix - 1
# df             <- data.frame(condition=c(rep(x='MET500',times=MET500.sample %>% length), rep(x='liver',times=GTex.liver.log2.read.count.matrix %>% ncol)))
# df$condition   <- factor(df$condition,levels = c('liver','MET500'))
# dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                          colData = df,
#                                          design = ~ condition)
# dds <- DESeq(dds)
# res1 <- results(dds,contrast = c('condition','MET500','liver')) %>% as.data.frame
# res1 <- res1[order(res1$pvalue),]
# res1 <- res1[complete.cases(res1),]
# 
# 
# 
# expressed.gene <- intersect(c(liver.expressed.gene,TCGA.expressed.gene) %>% unique,protein.coding.gene.id)
# expr.matrix    <- cbind(TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample],GTex.liver.log2.read.count.matrix[expressed.gene,])
# expr.matrix    <- 2^expr.matrix - 1
# df             <- data.frame(condition=c(rep(x='TCGA',times=TCGA.sample %>% length), rep(x='liver',times=GTex.liver.log2.read.count.matrix %>% ncol)))
# df$condition   <- factor(df$condition,levels = c('liver','TCGA'))
# dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                          colData = df,
#                                          design = ~ condition)
# dds <- estimateSizeFactors(dds)
# x   <-  counts(dds, normalized=TRUE)[APOH,]
# dds <- DESeq(dds)
# res2 <- results(dds,contrast = c('condition','liver','TCGA')) %>% as.data.frame
# res2 <- res2[order(res2$pvalue),]
# res2 <- res2[complete.cases(res2),]
# 
# 
# 
# expressed.gene <- intersect(c(MET500.expressed.gene,TCGA.expressed.gene) %>% unique,protein.coding.gene.id)
# expr.matrix    <- cbind(TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample],MET500.log2.read.count.matrix[expressed.gene,MET500.sample])
# expr.matrix    <- 2^expr.matrix - 1
# df             <- data.frame(condition=c(rep(x='TCGA',times=TCGA.sample %>% length), rep(x='MET500',times=MET500.sample %>% length)))
# df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
# dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                          colData = df,
#                                          design = ~ condition)
# dds <- estimateSizeFactors(dds)
# x   <-  counts(dds, normalized=TRUE)[APOH,]
# dds <- DESeq(dds)
# res3 <- results(dds,contrast = c('condition','MET500','TCGA')) %>% as.data.frame
# res3 <- res3[order(res3$pvalue),]
# res3 <- res3[complete.cases(res3),]
# 
# 
# 
# 
# 
# g <- intersect(TCGA.expressed.gene,MET500.expressed.gene)
# expressed.gene <- c(TCGA.expressed.gene,MET500.expressed.gene,liver.expressed.gene) %>% unique
# expr.matrix    <- cbind(TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample],MET500.log2.read.count.matrix[expressed.gene,MET500.sample],GTex.liver.log2.read.count.matrix[expressed.gene,])
# expr.matrix    <- 2^expr.matrix - 1
# df             <- data.frame(condition=c(rep(x='TCGA',times=TCGA.sample %>% length), rep(x='MET500',times=MET500.sample %>% length),rep(x='liver',times=GTex.liver.log2.read.count.matrix %>% ncol)   ))
# df$condition   <- factor(df$condition,levels = c('TCGA','MET500','liver'))
# dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                          colData = df,
#                                          design = ~ condition)
# dds <- estimateSizeFactors(dds)
# x   <-  counts(dds, normalized=TRUE)[APOH,]
# dds <- DESeq(dds)
# res4 <- results(dds,contrast = c('condition','MET500','TCGA')) %>% as.data.frame
# res4 <- res4[order(res4$pvalue),]
# res4 <- res4[complete.cases(res4),]
# 
# res5 <- results(dds,contrast = c('condition','liver','TCGA')) %>% as.data.frame
# res5 <- res5[order(res5$pvalue),]
# res5 <- res5[complete.cases(res5),]
# 
# 
# 
# ##############################
# g <- intersect(rownames(res1),rownames(res2))
# g <- intersect(rownames(res3),g)
# 
# 
# plot(x=res1[g,'log2FoldChange'],y=res3[g,'log2FoldChange'] - res2[g,'log2FoldChange'])
# lines(c(-20,20),c(-20,20))
# x=res1[g,'log2FoldChange']
# y=res3[g,'log2FoldChange'] - res2[g,'log2FoldChange']
# z <- y -x
# names(z) <- g
# 
# 
# 
# cut.off                        <- 0.05
# up.gene <- rownames(res)[res$log2FoldChange >  1  & res$padj < cut.off]
# dn.gene <- rownames(res)[res$log2FoldChange < -1  & res$padj < cut.off]
# 
# 
# up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)
# dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = dn.gene)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# TCGA.sample    <- pure.TCGA.breast.cancer.polyA.LumB.sample
# 
# 
# 
# 
# 
# expressed.gene <- intersect(c(breast.expressed.gene,liver.expressed.gene) %>% unique,protein.coding.gene.id)
# expr.matrix    <- cbind(GTex.breast.log2.read.count.matrix[expressed.gene,],GTex.liver.log2.read.count.matrix[expressed.gene,])
# expr.matrix    <- 2^expr.matrix - 1
# df             <- data.frame(condition=c(rep(x='breast',times=GTex.breast.log2.read.count.matrix %>% ncol), rep(x='liver',times=GTex.liver.log2.read.count.matrix %>% ncol)))
# df$condition   <- factor(df$condition,levels = c('breast','liver'))
# dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                          colData = df,
#                                          design = ~ condition)
# 
# dds <- DESeq(dds)
# res <- results(dds,contrast = c('condition','liver','breast')) %>% as.data.frame
# res <- res[order(res$pvalue),]
# res <- res[complete.cases(res),]
# cut.off                        <- 0.05
# de.res.liver.vs.breast.up.gene <- rownames(res)[res$log2FoldChange >  1  & res$padj < cut.off]
# de.res.liver.vs.breast.dn.gene <- rownames(res)[res$log2FoldChange < -1  & res$padj < cut.off]
# de.res.liver.vs.breast         <- res
