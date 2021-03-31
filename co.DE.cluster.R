library(org.Hs.eg.db)
require(dplyr)
require(bedr)
require(foreach)

gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id

load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')
DE.rs.list        <- list(BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, BRCA.LumB.DE.rs, PRAD.DE.rs,COAD.DE.rs,NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCA.Basal', 'BRCA.Her2', 'BRCA.LumB','PRAD', 'COAD', 'PNET', 'SINET')

load('client-side/output/hg19.gene.info.R.output/hg19.gene.info.RData')
distance <- 60000
get.co.up.DE.cluster <- function(DE.rs) {
    up.gene.symbol    <- mapping.df[DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,'symbol']
    up.gene.symbol    <- up.gene.symbol[is.na(up.gene.symbol) == FALSE]
    idx               <- match(up.gene.symbol,hg19.gene.info$genename)
    up.gene.coor.df   <- hg19.gene.info[idx,]     
    up.gene.coor.df   <- up.gene.coor.df[complete.cases(up.gene.coor.df),]

    bed.region        <- paste("chr",up.gene.coor.df$chrom,':',up.gene.coor.df$start,"-",up.gene.coor.df$end,sep='')
    names(bed.region) <- up.gene.coor.df$genename

    bed.sorted.region <- bedr.sort.region(bed.region)
    merged.region     <- bedr.merge.region(x=bed.sorted.region,distance=distance,check.chr = TRUE,check.sort = TRUE)

    
    co.up.clusters    <- setdiff(merged.region,bed.sorted.region)
    if(length(co.up.clusters) == 0){
        return(NULL)
    }
    region.size       <- bedr:::size.region(co.up.clusters)
    df                <- data.frame(co.up.clusters=co.up.clusters,region.size=region.size)
    df$co.up.clusters <- as.character(df$co.up.clusters)

    df$gene.name <- foreach(x= df$co.up.clusters %>% as.character(),.combine='c') %do% {
        f <- bed.sorted.region %in.region% x
        s <- bed.sorted.region[f]
        paste(names(bed.region)[match(s,bed.region)],collapse =":")
    }

    df$gene.number <- foreach(x= df$co.up.clusters %>% as.character(),.combine='c') %do% {
        f <- bed.sorted.region %in.region% x
        sum(f)
    }

    df <- df[order(df$gene.number,decreasing = TRUE),]
    df
}

get.co.dn.DE.cluster <- function(DE.rs) {
  dn.gene.symbol    <- mapping.df[DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene,'symbol']
  dn.gene.symbol    <- dn.gene.symbol[is.na(dn.gene.symbol) == FALSE]
  idx               <- match(dn.gene.symbol,hg19.gene.info$genename)
  dn.gene.coor.df   <- hg19.gene.info[idx,]     
  dn.gene.coor.df   <- dn.gene.coor.df[complete.cases(dn.gene.coor.df),]
  
  bed.region        <- paste("chr",dn.gene.coor.df$chrom,':',dn.gene.coor.df$start,"-",dn.gene.coor.df$end,sep='')
  names(bed.region) <- dn.gene.coor.df$genename
  
  bed.sorted.region <- bedr.sort.region(bed.region)
  merged.region     <- bedr.merge.region(x=bed.sorted.region,distance=distance,check.chr = TRUE,check.sort = TRUE)
  
  co.dn.clusters    <- setdiff(merged.region,bed.sorted.region)
  if(length(co.dn.clusters) == 0){
      return(NULL)
  }
  region.size       <- bedr:::size.region(co.dn.clusters)
  df                <- data.frame(co.dn.clusters=co.dn.clusters,region.size=region.size)
  df$co.dn.clusters <- as.character(df$co.dn.clusters)
  
  df$gene.name <- foreach(x= df$co.dn.clusters %>% as.character(),.combine='c') %do% {
    f <- bed.sorted.region %in.region% x
    s <- bed.sorted.region[f]
    paste(names(bed.region)[match(s,bed.region)],collapse =":")
  }
  
  df$gene.number <- foreach(x= df$co.dn.clusters %>% as.character(),.combine='c') %do% {
    f <- bed.sorted.region %in.region% x
    sum(f)
  }
  
  df <- df[order(df$gene.number,decreasing = TRUE),]
  df
}

co.up.DE.cluster.list <- lapply(DE.rs.list,get.co.up.DE.cluster)
co.dn.DE.cluster.list <- lapply(DE.rs.list,get.co.dn.DE.cluster)
names(co.up.DE.cluster.list) <- names(DE.rs.list)
names(co.dn.DE.cluster.list) <- names(DE.rs.list)


pooled.up.DE.cluster.df <- foreach(i = 1: length(co.up.DE.cluster.list),.combine='rbind') %do% {
    df             <-   co.up.DE.cluster.list[[i]]
    if( is.null(nrow(df)) == FALSE){
        df$cancer.type <- names(co.up.DE.cluster.list)[i]
        df
    }else{
        NULL  
    }
}

pooled.dn.DE.cluster.df <- foreach(i = 1: length(co.dn.DE.cluster.list),.combine='rbind') %do% {
  df             <-   co.dn.DE.cluster.list[[i]]
  if(is.null(nrow(df)) == FALSE){
    df$cancer.type <- names(co.dn.DE.cluster.list)[i]
    df
  }else{
    NULL  
  }
}

BRCA.LumB.fake.cluster.rs <- foreach(i = 1:1000, .combine='rbind') %do% {
    fake.DE.rs <- list()
    fake.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene <- sample(mapping.df$ensembl_gene_id,length(BRCA.LumB.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene))
    get.co.up.DE.cluster(fake.DE.rs)
}

BRCA.Basal.fake.cluster.rs <- foreach(i = 1:1000, .combine='rbind') %do% {
    fake.DE.rs <- list()
    fake.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene <- sample(mapping.df$ensembl_gene_id,length(BRCA.Basal.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene))
    get.co.up.DE.cluster(fake.DE.rs)
}

# PRAD.fake.cluster.rs <- foreach(i = 1:1000, .combine='rbind') %do% {
#     fake.DE.rs <- list()
#     fake.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene <- sample(mapping.df$ensembl_gene_id,length(PRAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene))
#     get.co.up.DE.cluster(fake.DE.rs)
# }


save(file = 'client-side/output/co.DE.cluster.R.output/co.DE.cluster.RData',list=c('pooled.up.DE.cluster.df','pooled.dn.DE.cluster.df','BRCA.LumB.fake.cluster.rs','BRCA.Basal.fake.cluster.rs'))


# chr19.p13.12.start <- 14000001
# chr19.p13.12.end   <- 16300000

start <- 13001942 #start of gene GCDH
end   <- 15560762 # end of  gene WIZ



flag              <- hg19.gene.info$chrom == '19' & hg19.gene.info$start >= start & hg19.gene.info$end <= end
chr19.p13.12.gene <- hg19.gene.info[flag,'genename'] %>% as.character()
flag              <- grepl(chr19.p13.12.gene,pattern = 'MIR') | grepl(chr19.p13.12.gene,pattern = 'LINC') | grepl(chr19.p13.12.gene,pattern = 'LOC') |grepl(chr19.p13.12.gene,pattern = 'SNORA')
chr19.p13.12.gene <- chr19.p13.12.gene[!flag]

basal.up.gene <- mapping.df[BRCA.Basal.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,'symbol'] %>% unique
basal.dn.gene <- mapping.df[BRCA.Basal.DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene,'symbol'] %>% unique

c.gene           <- intersect(basal.up.gene,chr19.p13.12.gene)
c.gene.df       <- hg19.gene.info[hg19.gene.info$genename %in% c.gene,]
c.gene.df       <- c.gene.df[order(c.gene.df$start),]
c.gene.df$chrom <- paste('chr',c.gene.df$chrom,sep='')
write.table(x=c.gene.df[,c('chrom','start','end','genename')],file='client-side/output/co.DE.cluster.R.output//basal.up.chr19.gene.bed',quote=FALSE,row.names=FALSE)



get.id <- function(x){
    l <- strsplit(x = x,split='\\.')  %>% unlist 
    l[1]
}
MET500.seg.data                         <- read.csv("~/Project/Cancer2CellLine/client-side/Data/CNV/cnv_v4.csv", stringsAsFactors=FALSE)
MET500.seg.data                         <- MET500.seg.data[,c('Pipeline_ID','Chr','Start','End','Log2_Coverage_Ratio')]
colnames(MET500.seg.data)               <- c('ID','chrom','loc.start','loc.end','seg.mean')
MET500.seg.data$ID                      <- sapply(MET500.seg.data$ID,get.id)

load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')

pure.MET.breast.cancer.Basal.patient <- MET500.sample.meta$MET500.id[match(pure.MET.breast.cancer.Basal.sample,MET500.sample.meta$Run)]
flag <-(MET500.seg.data$ID %in% pure.MET.breast.cancer.Basal.patient) & (MET500.seg.data$chrom == 19) &
       (MET500.seg.data$loc.end >= start) & (MET500.seg.data$loc.start <= end)
bed.data <- MET500.seg.data[flag,]  
bed.data$genename <- bed.data$ID
bed.data <- bed.data[order(bed.data$genename),]
write.table(x=bed.data[,c('chrom','loc.start','loc.end','genename')],file='client-side/output/co.DE.cluster.R.output//MET500.basal.cnv.call.bed',quote=FALSE,row.names=FALSE)






# ################################ analyze CNV data #############################
# load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')
# data      <- fread(input = 'client-side/Data/cBioPortal/brca_tcga_pan_can_atlas_2018/data_CNA.txt',header=TRUE) %>% as.data.frame
# col.id    <- match(pure.PRI.breast.cancer.Basal.sample,table = colnames(data))
# col.id    <- col.id[is.na(col.id) == FALSE]
# row.id    <- which(data$Hugo_Symbol == 'ZNF598') 
# 
# 
# data      <- fread(input = 'client-side/Data/cBioPortal/brca_tcga_pan_can_atlas_2018/data_log2CNA.txt',header=TRUE) %>% as.data.frame
# 

 
# 
# require(data.table)
# load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')
# cnv.data           <- fread('client-side/Data/cBioPortal/brca_tcga/data_linear_CNA.txt',header=TRUE) %>% as.data.frame
# rownames(cnv.data) <- cnv.data$Hugo_Symbol
# cnv.data.matrix    <- as.matrix(cnv.data[,c(-1,-2)])
# 
# 

# load('server-side/RData/Breast Invasive Carcinoma.RData')
# 
# # #NOTCH3.cnv <- TCGA.breast.cancer.cnv.matrix['NOTCH3',intersect(colnames(TCGA.breast.cancer.cnv.matrix),pure.PRI.breast.cancer.Basal.sample)]
# # #BRD4.cnv   <- TCGA.breast.cancer.cnv.matrix['BRD4',intersect(colnames(TCGA.breast.cancer.cnv.matrix),pure.PRI.breast.cancer.Basal.sample)]
# # RAD23A.cnv   <- TCGA.prostate.cancer.cnv.matrix['RECQL4',intersect(colnames(TCGA.prostate.cancer.cnv.matrix),pure.PRI.prostate.cancer.sample)]
# # 
# # RAD23A.cnv   <- MET500.prostate.cancer.cnv.matrix['RECQL4',intersect(colnames(MET500.prostate.cancer.cnv.matrix),MET500.sample.meta[pure.MET.prostate.cancer.sample,'MET500.id'])]
# load('client-side/output/organize.TCGA.and.MET500.breast.cancer.cnv.data.R.output/organize.TCGA.and.MET500.breast.cancer.cnv.data.RData')
# 
# cnv.data           <- fread('client-side/Data/cBioPortal/brca_tcga/data_linear_CNA.txt',header=TRUE) %>% as.data.frame
# rownames(cnv.data) <- cnv.data$Hugo_Symbol
# cnv.data.matrix    <- as.matrix(cnv.data[,c(-1,-2)])
# # MERTRN, good example of cnv- TME association
# cnv.vec.1            <- cnv.data.matrix[c('RAD23A'),intersect(colnames(cnv.data.matrix),c(pure.PRI.breast.cancer.Basal.sample))]
# cnv.vec.2            <- TCGA.breast.cancer.cnv.matrix['RAD23A',intersect(colnames(TCGA.breast.cancer.cnv.matrix),c(pure.PRI.breast.cancer.Basal.sample))]
# 
# # cnv.vec.1            <- cnv.data.matrix[c('LIN37'),intersect(colnames(cnv.data.matrix),c(pure.PRI.breast.cancer.LumB.sample,pure.PRI.breast.cancer.LumA.sample))]
# # cnv.vec.2            <- TCGA.breast.cancer.cnv.matrix['LIN37',intersect(colnames(TCGA.breast.cancer.cnv.matrix),c(pure.PRI.breast.cancer.LumB.sample,pure.PRI.breast.cancer.LumA.sample))]
# 
# 
# cnv.vec    <- cnv.vec.1
# TRE.sample <- names(cnv.vec)[cnv.vec >= 0.5]
# CON.sample <- names(cnv.vec)[cnv.vec < 0.5 ]
# 
# 
# # cnv.vec <- cnv.vec.1
# # TRE.sample <- names(cnv.vec)[cnv.vec >= 0.5]
# # CON.sample <- names(cnv.vec)[cnv.vec < 0.5 ]
# # 
# 
# 
# 
# rs <- perform.DE.analysis.between.TRE.and.CON(CON.log2.tpm.matrix         = log2.tpm.matrix[,CON.sample],
#                                               CON.log2.read.count.matrix  = log2.read.count.matrix[,CON.sample],
#                                               TRE.log2.tpm.matrix         = log2.tpm.matrix[,TRE.sample],
#                                               TRE.log2.read.count.matrix  = log2.read.count.matrix[,TRE.sample]
# 
# )
# up.gene.1            <- rownames(rs)[rs$log2FoldChange > 0.5  & rs$padj < 0.05]
# dn.gene.1            <- rownames(rs)[rs$log2FoldChange < -0.5 & rs$padj < 0.05]
# 
# 
# 
# T.cell.marker.gene <- c('ENSG00000153563','ENSG00000172116','ENSG00000145649','ENSG00000100453','ENSG00000180644')
# s <- names(cnv.vec)
# 
# m <- 2^log2.tpm.matrix[T.cell.marker.gene,s] - 1
# 
# gm <- function(x) {
#     exp(mean(log(x)))
# }
# 
# CD8.T.cell.score <- apply(m,2,gm)
# CD8.T.cell.score <- CD8.T.cell.score[is.na(CD8.T.cell.score) == FALSE]
# 
# tmp <- intersect(names(CD8.T.cell.score),names(cnv.vec))
# 
# plot(x=cnv.vec[tmp],y=CD8.T.cell.score[tmp] %>% log2)
# 
# 
# TRE.sample <- intersect(TRE.sample,tmp)
# CON.sample <- intersect(CON.sample,tmp)
# 
# 
# boxplot(CD8.T.cell.score[CON.sample] %>% log,CD8.T.cell.score[TRE.sample] %>% log)
# 
# t.test(CD8.T.cell.score[CON.sample] %>% log,CD8.T.cell.score[TRE.sample] %>% log)
# 
# ks.test(CD8.T.cell.score[CON.sample] %>% log,CD8.T.cell.score[TRE.sample] %>% log)
# 
# wilcox.test(CD8.T.cell.score[CON.sample] %>% log,CD8.T.cell.score[TRE.sample] %>% log)
# 
# 
# load('~/Project/InSilicoCRISPR/client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData')
# 
# purity <- TCGA.CD8.T.cell.level.df[names(CD8.T.cell.score),'tumor.purity']
# names(purity) <- names(CD8.T.cell.score)
# boxplot(purity[CON.sample],purity[TRE.sample])
# 
# df <- data.frame(purity=purity,score = CD8.T.cell.score %>% log)
# df <- df[complete.cases(df),]
# loess.fit.rs     <- loess(data = df,formula= score ~ purity,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# adjusted.score   <- loess.fit.rs$residuals
# 
# plot(x=purity[names(adjusted.score)],y=adjusted.score)
# 
# 
# 
# plot(x=cnv.vec[names(adjusted.score)],y=(adjusted.score) )
# 
# TRE.sample <- names(cnv.vec)[cnv.vec >= 1]
# CON.sample <- names(cnv.vec)[cnv.vec < 0.3 & cnv.vec > -0.3 ]
# boxplot(adjusted.score[CON.sample],adjusted.score[TRE.sample])
# ks.test(adjusted.score[CON.sample],adjusted.score[TRE.sample])
# wilcox.test(adjusted.score[CON.sample],adjusted.score[TRE.sample])
# t.test(adjusted.score[CON.sample],adjusted.score[TRE.sample])
# 
# 
# 
# 
# 
# 
# 
# 
# 
# H.sample <- intersect(names(cnv.vec)[cnv.vec >= 0.5],tmp)
# 
# y <- sort(CD8.T.cell.score[tmp] %>% log2)
# 
# plot(1:length(y),y=ifelse(names(y) %in% H.sample,1,0),col=ifelse(names(y) %in% H.sample,'red','black'),pch=19)
# 
# require(fgsea)
# rs <- fgsea(pathways = list(H.sample=H.sample),stats = (y),nperm = 10000)
# plotEnrichment(pathway = H.sample,rank(y))
