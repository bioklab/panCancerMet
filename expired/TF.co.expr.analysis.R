source('client-side/code/DEBoost.R')


load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')


DE.rs.list        <- list(BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, BRCA.LumB.DE.rs, COAD.DE.rs, PRAD.DE.rs, NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCA.Basal', 'BRCA.Her2', 'BRCA.LumB', 'COAD','PRAD', 'NET.PAAD', 'NET.SI')


gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id


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

TF.mapping.df <- mapping.df[mapping.df$symbol %in% TF.list,]
TF.ensmeble.id <- TF.mapping.df$ensembl_gene_id %>% as.character()

get.up.TF.co.expr.gene <- function(DE.rs) {
   up.DE.gene  <-  DE.rs$tumor.intrinsic.DE.gene.rs$up.gene  
   up.DE.TF    <-  intersect(up.DE.gene,TF.ensmeble.id)
   up.DE.gene  <- setdiff(up.DE.gene,up.DE.TF)
   cor.matrix  <- cor(PRI.log2.tpm.matrix[up.DE.TF,] %>% t,  PRI.log2.tpm.matrix[up.DE.gene,] %>% t,method='spearman')
   idx         <- apply(cor.matrix, 2, function(x) ifelse( sum(x >= 0.5) > 0, TRUE,FALSE) ) %>% which
   tmp         <- cor.matrix[,idx]
   co.expr.gene.vec <- apply(cor.matrix,1,function(x) sum(x >= 0.5) ) %>% sort(decreasing = TRUE)
   tmp         <- tmp[names(co.expr.gene.vec),]
   flag <- which(tmp >= 0.5,arr.ind = 2)
   tmp[flag] <- 1
   flag <- which(tmp < 0.5,arr.ind = 2)
   tmp[flag] <- 0
   tmp
}


get.dn.TF.co.expr.gene <- function(DE.rs) {
    up.DE.gene  <-  DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene  
    up.DE.TF    <-  intersect(up.DE.gene,TF.ensmeble.id)
    up.DE.gene  <- setdiff(up.DE.gene,up.DE.TF)
    cor.matrix  <- cor(PRI.log2.tpm.matrix[up.DE.TF,] %>% t,  PRI.log2.tpm.matrix[up.DE.gene,] %>% t,method='spearman')
    idx         <- apply(cor.matrix, 2, function(x) ifelse( sum(x >= 0.5) > 0, TRUE,FALSE) ) %>% which
    tmp         <- cor.matrix[,idx]
    co.expr.gene.vec <- apply(cor.matrix,1,function(x) sum(x >= 0.5) ) %>% sort(decreasing = TRUE)
    tmp         <- tmp[names(co.expr.gene.vec),]
    flag <- which(tmp >= 0.5,arr.ind = 2)
    tmp[flag] <- 1
    flag <- which(tmp < 0.5,arr.ind = 2)
    tmp[flag] <- 0
    tmp

}


load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')

PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.Basal.sample]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]
m.BRCA.Basal               <- get.up.TF.co.expr.gene(DE.rs.list$BRCA.Basal)
n.BRCA.Basal               <- get.dn.TF.co.expr.gene(DE.rs.list$BRCA.Basal)




load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.Her2.sample]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,pure.PRI.breast.cancer.Her2.sample]
m <- get.up.TF.co.expr.gene(DE.rs.list$BRCA.Her2)
n <- get.dn.TF.co.expr.gene(DE.rs.list$BRCA.Her2)





load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.LumB.sample]
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.LumB.sample]
m <- get.up.TF.co.expr.gene(DE.rs.list$BRCA.LumB)
n <- get.dn.TF.co.expr.gene(DE.rs.list$BRCA.LumB)


load('client-side/output/Select.pure.sample.colorectal.cancer.R.output/Select.pure.sample.colorectal.cancer.RData')
load('server-side/RData/COLORECTAL_SRP029880.RData')
PRI.log2.read.count.matrix <- COLORECTAL_SRP029880_log2.read.count.matrix[,pure.PRI.colorectal.cancer.sample]
PRI.log2.tpm.matrix        <- COLORECTAL_SRP029880_log2.tpm.matrix[,pure.PRI.colorectal.cancer.sample]
co.expr.rs                 <- c(get.up.TF.co.expr.gene(DE.rs.list$COAD) ,get.dn.TF.co.expr.gene(DE.rs.list$COAD))


load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')
load('server-side/RData//Prostate Adenocarcinoma.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.prostate.cancer.sample]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]
m.PRAD                     <- get.up.TF.co.expr.gene(DE.rs.list$PRAD)
n.PRAD                     <- get.dn.TF.co.expr.gene(DE.rs.list$PRAD)


load('client-side/output/Select.pure.sample.NET.pancreatic.cancer.R.output/Select.pure.sample.NET.pancreatic.cancer.RData')
load('server-side/RData/GEP.NET.RData')
PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
PRI.log2.tpm.matrix        <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
m <- get.up.TF.co.expr.gene(DE.rs.list$NET.PAAD)
n <- get.dn.TF.co.expr.gene(DE.rs.list$NET.PAAD)


load('client-side/output/Select.pure.sample.NET.si.cancer.R.output/Select.pure.sample.NET.si.cancer.RData')
load('server-side/RData/GEP.NET.RData')
PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.si.cancer.sample]
PRI.log2.tpm.matrix       <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.si.cancer.sample]
m <- get.up.TF.co.expr.gene(DE.rs.list$NET.SI)
n <- get.dn.TF.co.expr.gene(DE.rs.list$NET.SI)








gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id


up.gene                  <- lapply(TWIST1.DE.rs.list,function(x) intersect(x$up.gene, mapping.df$ensembl_gene_id)    ) %>% unlist %>% unique
up.gene.matrix           <- matrix(0,nrow = length(up.gene),ncol = length(TWIST1.DE.rs.list))
rownames(up.gene.matrix) <- up.gene
colnames(up.gene.matrix) <- names(TWIST1.DE.rs.list)
for(i in 1:length(TWIST1.DE.rs.list) )  {
  gene   <- intersect(TWIST1.DE.rs.list[[i]]$up.gene, mapping.df$ensembl_gene_id)
  cancer <- names(TWIST1.DE.rs.list)[i]
  up.gene.matrix[gene,cancer] <- 1
}
up.gene.freq <- apply(up.gene.matrix,1,sum)
TWIST1.up.gene.freq <- up.gene.freq

dn.gene                  <- lapply(TWIST1.DE.rs.list,function(x) intersect(x$dn.gene, mapping.df$ensembl_gene_id)) %>% unlist %>% unique
dn.gene.matrix           <- matrix(0,nrow = length(dn.gene),ncol = length(TWIST1.DE.rs.list))
rownames(dn.gene.matrix) <- dn.gene
colnames(dn.gene.matrix) <- names(TWIST1.DE.rs.list)
for(i in 1:length(TWIST1.DE.rs.list) )  {
  gene   <- intersect(TWIST1.DE.rs.list[[i]]$dn.gene,mapping.df$ensembl_gene_id)
  cancer <- names(TWIST1.DE.rs.list)[i]
  dn.gene.matrix[gene,cancer] <- 1
}
dn.gene.freq <- apply(dn.gene.matrix,1,sum)
TWIST1.dn.gene.freq <- dn.gene.freq

TWIST1.dn.signature <- names(TWIST1.dn.gene.freq)[TWIST1.dn.gene.freq >= 4]
TWIST1.up.signature <- names(TWIST1.dn.gene.freq)[TWIST1.up.gene.freq >= 4]

