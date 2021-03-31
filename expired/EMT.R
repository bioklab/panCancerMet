require(dplyr)
require(stringr)

file.name <- "client-side/Data/JASPAR.txt"
conn      <- file(file.name,open="r")
line      <- readLines(conn)
get.TF.name <- function(x) {
  if(grepl(x = x,pattern='>')){
    l <- strsplit(x=x,split='\t')  %>% unlist
    toupper(l[2])
  }else{
    NA
  }
}

TF.list <- sapply(line,get.TF.name)
TF.list <- TF.list[is.na(TF.list) == FALSE]
names(TF.list) <- NULL
JASPAR.TF.list <- TF.list
close(conn)


file.name <- "client-side/Data/CISTROME.factor.line.txt"
conn      <- file(file.name,open="r")
line      <- readLines(conn)
tmp       <- str_extract_all(pattern="id=\"[:alnum:]+\"",string=line[1],simplify = TRUE)
get.TF.name <- function(x) {
  x <- str_remove_all(string = x,pattern="\"")
  x <- str_remove_all(string = x,pattern="id=")
  x
}
TF.list <- sapply(tmp,get.TF.name)
TF.list <- TF.list[is.na(TF.list) == FALSE]
names(TF.list) <- NULL
CISTROME.TF.list <- TF.list
close(conn)

TF.list <- c(CISTROME.TF.list,JASPAR.TF.list) %>% unique


SNAI1 <- 'ENSG00000124216'
SNAI2 <- 'ENSG00000019549'
ZEB1  <- 'ENSG00000148516'
ZEB2  <- 'ENSG00000169554'
TWIST <- 'ENSG00000122691'
CDH1  <- 'ENSG00000039068'
EMT.TF <- c(SNAI1,SNAI2,ZEB1,ZEB2,TWIST,CDH1)

basal.up.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
basal.dn.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
intersect(basal.up.gene.symbol,TF.list)
intersect(basal.dn.gene.symbol,TF.list)

lumb.up.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.up.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
lumb.dn.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
intersect(lumb.up.gene.symbol,TF.list)
intersect(lumb.dn.gene.symbol,TF.list)

her2.up.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.up.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
her2.dn.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
intersect(her2.up.gene.symbol,TF.list)
intersect(her2.dn.gene.symbol,TF.list)


x1 <- intersect(basal.dn.gene.symbol,TF.list)
x2 <- intersect(lumb.dn.gene.symbol, TF.list)
x3 <- intersect(her2.dn.gene.symbol, TF.list)
tmp.dn.TF <- table(c(x1,x2,x3)) %>% as.data.frame



x1 <- intersect(basal.up.gene.symbol,TF.list)
x2 <- intersect(lumb.up.gene.symbol, TF.list)
x3 <- intersect(her2.up.gene.symbol, TF.list)
tmp.up.TF <- table(c(x1,x2,x3)) %>% as.data.frame



load('server-side/RData/Breast Invasive Carcinoma.RData')
load('client-side/output/TCGA.breast.cancer.meta.R.output/TCGA.breast.cancer.meta.RData')
get.outlier <- function(x) {
    q <- quantile(x)  
    upper <- q[4] + 1.5 * IQR(x)
    lower <- q[1] - 1.5 * IQR(x)
    list(upper=names(x)[x > upper],lower=names(x)[x < lower])
}
PRRX1 <- 'ENSG00000116132'

basal.dn.gene.id <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$ENSEMBL
cor.matrix       <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.Basal.sample] ,log2.fpkm.matrix[,pure.TCGA.breast.cancer.polyA.Basal.sample] %>% t,method='spearman')
cor.vec          <- c(cor.matrix)
names(cor.vec)   <- rownames(log2.fpkm.matrix)
cor.vec          <- sort(cor.vec,decreasing = TRUE)
cor.vec          <- cor.vec[is.na(cor.vec) == FALSE]
upper.gene       <- get.outlier(cor.vec)$upper





lumb.dn.gene.id <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$ENSEMBL
cor.matrix       <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.LumB.sample] ,log2.fpkm.matrix[,pure.TCGA.breast.cancer.polyA.LumB.sample] %>% t,method='spearman')
cor.vec          <- c(cor.matrix)
names(cor.vec)   <- rownames(log2.fpkm.matrix)
cor.vec          <- sort(cor.vec,decreasing = TRUE)
cor.vec          <- cor.vec[is.na(cor.vec) == FALSE]
upper.gene       <- get.outlier(cor.vec)$upper
wilcox.test(cor.vec,cor.vec[lumb.dn.gene.id])




her2.dn.gene.id  <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$ENSEMBL
cor.matrix       <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.Her2.sample] ,log2.fpkm.matrix[,pure.TCGA.breast.cancer.polyA.Her2.sample] %>% t,method='spearman')
cor.vec          <- c(cor.matrix)
names(cor.vec)   <- rownames(log2.fpkm.matrix)
cor.vec          <- sort(cor.vec,decreasing = TRUE)
cor.vec          <- cor.vec[is.na(cor.vec) == FALSE]
upper.gene       <- get.outlier(cor.vec)$upper
wilcox.test(cor.vec,cor.vec[her2.dn.gene.id])
boxplot(cor.vec,cor.vec[her2.dn.gene.id])





