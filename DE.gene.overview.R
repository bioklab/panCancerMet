require(clusterProfiler)
library(org.Hs.eg.db)
require(dplyr)
require(pheatmap)

gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id


load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')
load('client-side/output/DE.prostate.cancer.SRP253428.R.output/DE.prostate.cancer.SRP253428.RData')

DE.rs.list        <- list(BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, BRCA.LumB.DE.rs, PRAD.DE.rs, COAD.DE.rs, NET.PAAD.DE.rs, NET.SI.DE.rs)
names(DE.rs.list) <- c(   'BRCA.Basal',     'BRCA.Her2',     'BRCA.LumB',     'PRAD',     'COAD',     'PNET',     'SINET')


perform.GO.analysis <- function(DE.rs) {
    up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = intersect(DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,mapping.df$ensembl_gene_id))
    up.gene.annotation.df  <- up.gene.annotation.df[complete.cases(up.gene.annotation.df),]
    GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL %>% unique ,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
    GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 , ]
    GO.rs.1.up             <- GO.rs.1.up[order(GO.rs.1.up$pvalue),]
    
    dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = intersect(DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene,mapping.df$ensembl_gene_id))
    dn.gene.annotation.df  <- dn.gene.annotation.df[complete.cases(dn.gene.annotation.df),]
    GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL %>% unique,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
    GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 , ]
    GO.rs.1.dn             <- GO.rs.1.dn[order(GO.rs.1.dn$pvalue),]
    
    list(up = GO.rs.1.up,dn = GO.rs.1.dn)
}


GO.rs.list <- foreach(DE.rs = DE.rs.list) %do% {
    perform.GO.analysis(DE.rs)
}
names(GO.rs.list) <- names(DE.rs.list)









save(file='client-side/output/DE.gene.overview.R.output/DE.gene.overview.RData',list=c('GO.rs.list'))

# rs1 <- perform.GO.analysis(BRCA.LumB.DE.rs)
# rs2 <- perform.GO.analysis(BRCA.Basal.DE.rs)
# rs3 <- perform.GO.analysis(BRCA.Her2.DE.rs)
# rs4 <- perform.GO.analysis(PRAD.DE.rs)
# rs5 <- perform.GO.analysis(COAD.DE.rs)
# rs6 <- perform.GO.analysis(NET.PAAD.DE.rs)
# rs7 <- perform.GO.analysis(NET.SI.DE.rs)


compute.jaccard.index  <- function(x,y){
  ( intersect(x,y) %>% length )/( unique(c(x,y))   %>% length )
  
}

jacard.score.matrix <- matrix(data=0,nrow=length(DE.rs.list),ncol=length(DE.rs.list))
for(i in 1:length(DE.rs.list)) {
   for(j in 1:length(DE.rs.list)){
     if(i < j ){
       jacard.score.matrix[i,j] <- compute.jaccard.index(DE.rs.list[[i]]$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs.list[[j]]$tumor.intrinsic.DE.gene.rs$up.gene)
     }
     if( i == j){
       jacard.score.matrix[i,j] <- NA
     }
     if(i > j) {
       jacard.score.matrix[i,j] <- compute.jaccard.index(DE.rs.list[[i]]$tumor.intrinsic.DE.gene.rs$dn.gene,
                                                         DE.rs.list[[j]]$tumor.intrinsic.DE.gene.rs$dn.gene
                                                         )
       
     }
   }  
  
  
}
rownames(jacard.score.matrix) <- names(DE.rs.list)
colnames(jacard.score.matrix) <- names(DE.rs.list)
pheatmap(jacard.score.matrix[1:7,1:7],cluster_rows = F,cluster_cols = F,na_col = 'black',display_numbers = TRUE,number_color = 'black',fontsize_number = 25,number_format = '%.3f')


lt.vec <- c()
ut.vec <- c()
for(i in 1:(length(DE.rs.list) - 1) ){
    for(j in (i + 1): length(DE.rs.list)) { 
        lt.vec <- c(lt.vec, jacard.score.matrix[j,i])
        ut.vec <- c(ut.vec, jacard.score.matrix[i,j])
    }
}
wilcox.test(lt.vec,ut.vec,paired=TRUE)

s <- dn.gene.matrix[,c('BRCA.Basal','BRCA.LumB','BRCA.Her2','PRAD')]
h <- apply(s,1,sum)
batch.dn.gene <- names(h)[h == 4]

s <- up.gene.matrix[,c('BRCA.Basal','BRCA.LumB','BRCA.Her2','PRAD')]
h <- apply(s,1,sum)
batch.up.gene <- names(h)[h == 4]



DE.rs.list.copy <- DE.rs.list

DE.rs.list$BRCA.Basal$tumor.intrinsic.DE.gene.rs$dn.gene <- setdiff(DE.rs.list$BRCA.Basal$tumor.intrinsic.DE.gene.rs$dn.gene,bch.gene)
DE.rs.list$BRCA.Her2$tumor.intrinsic.DE.gene.rs$dn.gene  <- setdiff(DE.rs.list$BRCA.Her2$tumor.intrinsic.DE.gene.rs$dn.gene,bch.gene)
DE.rs.list$BRCA.LumB$tumor.intrinsic.DE.gene.rs$dn.gene  <- setdiff(DE.rs.list$BRCA.LumB$tumor.intrinsic.DE.gene.rs$dn.gene,bch.gene)
DE.rs.list$PRAD$tumor.intrinsic.DE.gene.rs$dn.gene       <- setdiff(DE.rs.list$PRAD$tumor.intrinsic.DE.gene.rs$dn.gene,bch.gene)

DE.rs.list$BRCA.Basal$tumor.intrinsic.DE.gene.rs$up.gene <- setdiff(DE.rs.list$BRCA.Basal$tumor.intrinsic.DE.gene.rs$up.gene,bch.gene)
DE.rs.list$BRCA.Her2$tumor.intrinsic.DE.gene.rs$up.gene  <- setdiff(DE.rs.list$BRCA.Her2$tumor.intrinsic.DE.gene.rs$up.gene,bch.gene)
DE.rs.list$BRCA.LumB$tumor.intrinsic.DE.gene.rs$up.gene  <- setdiff(DE.rs.list$BRCA.LumB$tumor.intrinsic.DE.gene.rs$up.gene,bch.gene)
DE.rs.list$PRAD$tumor.intrinsic.DE.gene.rs$up.gene       <- setdiff(DE.rs.list$PRAD$tumor.intrinsic.DE.gene.rs$up.gene,bch.gene)





jacard.score.matrix <- matrix(data=0,nrow=length(DE.rs.list),ncol=length(DE.rs.list))
for(i in 1:length(DE.rs.list)) {
  for(j in 1:length(DE.rs.list)){
    if(i < j ){
      jacard.score.matrix[i,j] <- compute.jaccard.index(DE.rs.list[[i]]$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs.list[[j]]$tumor.intrinsic.DE.gene.rs$up.gene)
    }
    if( i == j){
      jacard.score.matrix[i,j] <- NA
    }
    if(i > j) {
      jacard.score.matrix[i,j] <- compute.jaccard.index(DE.rs.list[[i]]$tumor.intrinsic.DE.gene.rs$dn.gene,
                                                        DE.rs.list[[j]]$tumor.intrinsic.DE.gene.rs$dn.gene
      )
      
    }
  }  
  
  
}
rownames(jacard.score.matrix) <- names(DE.rs.list)
colnames(jacard.score.matrix) <- names(DE.rs.list)
pheatmap(jacard.score.matrix,cluster_rows = F,cluster_cols = F,na_col = 'black',display_numbers = TRUE,number_color = 'black',fontsize_number = 25,number_format = '%.3f')



DE.rs.list <- DE.rs.list.copy

# g.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = 　c.dn.gene)
# g.df <- g.df[complete.cases(g.df),]
# GO.rs             <- enrichGO(gene=g.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
# GO.rs             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 & GO.rs.1.up$pvalue < 0.01, ]
# 



up.gene                  <- lapply(DE.rs.list,function(x) x$tumor.intrinsic.DE.gene.rs$up.gene ) %>% unlist %>% unique
up.gene.matrix           <- matrix(0,nrow = length(up.gene),ncol = length(DE.rs.list))
rownames(up.gene.matrix) <- up.gene
colnames(up.gene.matrix) <- names(DE.rs.list)
for(i in 1:length(DE.rs.list) )  {
    gene   <- DE.rs.list[[i]]$tumor.intrinsic.DE.gene.rs$up.gene
    cancer <- names(DE.rs.list)[i]
    up.gene.matrix[gene,cancer] <- 1
}
up.gene.freq <- apply(up.gene.matrix,1,sum)
up.gene.freq <- sort(up.gene.freq,decreasing = TRUE)
h.up.gene <- names(up.gene.freq)[up.gene.freq >= 2]
# pheatmap(up.gene.matrix[h.up.gene,])
# 
# 
# 
dn.gene                  <- lapply(DE.rs.list,function(x) x$tumor.intrinsic.DE.gene.rs$dn.gene ) %>% unlist %>% unique
dn.gene.matrix           <- matrix(0,nrow = length(dn.gene),ncol = length(DE.rs.list))
rownames(dn.gene.matrix) <- dn.gene
colnames(dn.gene.matrix) <- names(DE.rs.list)
for(i in 1:length(DE.rs.list) )  {
  gene   <- DE.rs.list[[i]]$tumor.intrinsic.DE.gene.rs$dn.gene
  cancer <- names(DE.rs.list)[i]
  dn.gene.matrix[gene,cancer] <- 1
}
dn.gene.freq <- apply(dn.gene.matrix,1,sum)
dn.gene.freq <- sort(dn.gene.freq,decreasing = TRUE)
h.dn.gene <- names(dn.gene.freq)[dn.gene.freq >= 4]






median.tpm.matrix <- read.gct(file='client-side/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
get.gene.id <- function(x) {
  strsplit(x=x,split='\\.') %>% unlist %>% head(1)
}
rownames(median.tpm.matrix) <- sapply(rownames(median.tpm.matrix),get.gene.id)

is.immune.gene <- function(x) {
  y    <- log2(x+1)    
  tmp  <- quantile(y)
  upper <- tmp[4] + 1.5* IQR(y)
  outlier.tissue <- names(y)[ y > upper]  
  immune.tissue    <- c('Spleen','Whole.Blood')
  if(sum(outlier.tissue %in% immune.tissue) ==2 & ('Liver' %in% outlier.tissue == FALSE) ){ # well, let us exclude liver specific genes
    return(TRUE)
  }else{
    return(FALSE)  
  }
}
flag               <- apply(median.tpm.matrix,1,is.immune.gene)
immune.gene        <- rownames(median.tpm.matrix)[flag]
names(immune.gene) <- NULL
write.csv(x=immune.gene,file='client-side/output/analyze.DE.gene.R.output/immune.gene.list.csv',quote=FALSE)


up.de.immune.gene.list <- foreach(i = 1:length(DE.rs.list)) %do% {
    obj <- DE.rs.list[[i]]$deseq2.M.vs.P.res
    up.de.gene <- rownames(obj)[obj$log2FoldChange > 1 & obj$padj < 0.05 ]
    intersect(up.de.gene,immune.gene)
}
names(up.de.immune.gene.list) <- names(DE.rs.list)
dn.de.immune.gene.list <- foreach(i = 1:length(DE.rs.list)) %do% {
    obj <- DE.rs.list[[i]]$deseq2.M.vs.P.res
    dn.de.gene <- rownames(obj)[obj$log2FoldChange < -1 & obj$padj < 0.05 ]
    intersect(dn.de.gene,immune.gene)
}
names(dn.de.immune.gene.list) <- names(DE.rs.list)




# pheatmap(dn.gene.matrix[h.dn.gene,])
# 
# pheatmap(dn.gene.matrix)
# 
# 
# CDH1 <- 'ENSG00000039068'
# CDH2 <- 'ENSG00000170558'
# FN1  <- 'ENSG00000115414'
# VIM  <- 'ENSG00000026025'




########### Analyze TF ################

# require(dplyr)
# require(stringr)
# 
# file.name <- "client-side/Data/JASPAR.txt"
# conn      <- file(file.name,open="r")
# line      <- readLines(conn)
# get.TF.name <- function(x) {
#   if(grepl(x = x,pattern='>')){
#     l <- strsplit(x=x,split='\t')  %>% unlist
#     toupper(l[2])
#   }else{
#     NA
#   }
# }
# 
# TF.list <- sapply(line,get.TF.name)
# TF.list <- TF.list[is.na(TF.list) == FALSE]
# names(TF.list) <- NULL
# JASPAR.TF.list <- TF.list
# close(conn)
# 
# 
# file.name <- "client-side/Data/CISTROME.factor.line.txt"
# conn      <- file(file.name,open="r")
# line      <- readLines(conn)
# tmp       <- str_extract_all(pattern="id=\"[:alnum:]+\"",string=line[1],simplify = TRUE)
# get.TF.name <- function(x) {
#   x <- str_remove_all(string = x,pattern="\"")
#   x <- str_remove_all(string = x,pattern="id=")
#   x
# }
# TF.list <- sapply(tmp,get.TF.name)
# TF.list <- TF.list[is.na(TF.list) == FALSE]
# names(TF.list) <- NULL
# CISTROME.TF.list <- TF.list
# close(conn)
# 
# TF.list <- c(CISTROME.TF.list,JASPAR.TF.list) %>% unique
# 
# up.TF.list <- lapply(DE.rs.list,function(x) intersect(mapping.df[x$tumor.intrinsic.DE.gene.rs$up.gene,'symbol'],TF.list)    )
# dn.TF.list <- lapply(DE.rs.list,function(x) intersect(mapping.df[x$tumor.intrinsic.DE.gene.rs$dn.gene,'symbol'],TF.list)    )
# 
# 
# up.TF.list.ensemble.id <-  mapping.df$ensembl_gene_id[match(x=unlist(up.TF.list) %>% unique,table = mapping.df$symbol)]
# dn.TF.list.ensemble.id <-  mapping.df$ensembl_gene_id[match(x=unlist(dn.TF.list) %>% unique,table = mapping.df$symbol)]
# 
# 
# 
# all.TF.gene.ensemble.id <- mapping.df$ensembl_gene_id[match(x=(TF.list) %>% unique,table = mapping.df$symbol)]
# 
# 
# 
# 
# 
# all.DE.gene <- c(BRCA.Basal.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene)
# DE.TF       <- intersect(all.TF.gene.ensemble.id,all.DE.gene)
# DE.TF       <- intersect(DE.TF, rownames(PRI.log2.tpm.matrix))
# 
# cor.matrix <- cor(PRI.log2.tpm.matrix[DE.TF,] %>% t,PRI.log2.tpm.matrix[setdiff(all.DE.gene,DE.TF),] %>% t,method='spearman')
# name.vec 
# 
# l <- foreach(i = 1:nrow(cor.matrix),.combine='rbind') %do% {
#     vec <- cor.matrix[i,]  
#     w   <- which(vec > 0.5)
#     if(sum(w) > 0){
#         data.frame(size = length(w),m.cor=median(cor.matrix[i,w]))
#     }else{
#         data.frame(size = -1,m.cor=-1)
# 
#     }
#   
#   
# }
# rownames(l) <- rownames(cor.matrix)
# 
# cc <- c('ENSG00000104885',names(w))
# 
# dd <- setdiff(all.DE.gene,cc)
# 
# 
# all.cor.matrix <- cor(PRI.log2.tpm.matrix[all.DE.gene,] %>% t,method='spearman')
# pheatmap(all.cor.matrix[c(cc,dd),c(cc,dd)] %>% abs ,cluster_rows = F,cluster_cols = F)
# 
# 
# 
# 
# 
# 
# 
# tmp <- unlist(l)
# 
# n <- sapply(l,function(x) x %>% length)
# 
# cor.tmp <- cor(PRI.log2.tpm.matrix[tmp,] %>% t %>% abs,method='spearman')
# 'ENSG00000104885: DOT1l'
# mapping.df[l$ENSG00000141867,]
# 




# 
# 
# c.up.gene        <- names(up.gene.freq)[up.gene.freq >=4]
# c.up.gene.symbol <- mapping.df[c.up.gene,'symbol']
# c.up.TF          <- intersect(c.up.gene.symbol,TF.list)
# 
# 
# 
# 
# c.dn.gene        <- names(dn.gene.freq)[dn.gene.freq >=4]
# c.dn.gene.symbol <- mapping.df[c.dn.gene,'symbol']
# c.dn.TF          <- intersect(c.dn.gene.symbol,TF.list)
# 
# 
# 
# 
# 
# 
# cosine_sim <- function(a, b) crossprod(a,b)/sqrt(crossprod(a)*crossprod(b))
# 
# 
# 
# 
# 
# 
# 
# PRRX1 <- "ENSG00000116132" # now, focus on PRRX1
# EGR1  <- "ENSG00000120738" # now, focus on EGR1
# SPARCL1 <- "ENSG00000152583"
# 
# 
# 
# PRRX1.expr   <- PRI.log2.fpkm.matrix[PRRX1,] 
# q            <- quantile(PRRX1.expr)
# h.sample     <- names(PRRX1.expr)[PRRX1.expr >= q['75%']]
# l.sample     <- names(PRRX1.expr)[PRRX1.expr <= q['25%']]
# 
# 
# rs <- perform.DE.analysis.between.TRE.and.CON(CON.log2.fpkm.matrix = PRI.log2.fpkm.matrix[,h.sample],
#                                               TRE.log2.fpkm.matrix = PRI.log2.fpkm.matrix[,l.sample],
#                                               CON.log2.read.count.matrix = PRI.log2.read.count.matrix[,h.sample],
#                                               TRE.log2.read.count.matrix = PRI.log2.read.count.matrix[,l.sample]
#                                               )
# 
# up.gene <- rownames(rs)[rs$log2FoldChange >  1 & rs$padj < 0.05]
# dn.gene <- rownames(rs)[rs$log2FoldChange < -1 & rs$padj < 0.05]
# 
# 
# w.test.p.value <- foreach(g=c(up.gene,dn.gene),.combine='c') %do% {
#     wilcox.test(PRI.log2.fpkm.matrix[g,h.sample],PRI.log2.fpkm.matrix[g,l.sample])$p.value  
# }
# names(w.test.p.value) <- c(up.gene,dn.gene)
# 
# de.gene <- names(w.test.p.value)[w.test.p.value < 0.05]
# 
# # m1 <- apply(PRI.log2.fpkm.matrix[de.gene,h.sample],1,median)
# # m2 <- apply(PRI.log2.fpkm.matrix[de.gene,l.sample],1,median)
# # 
# # max.m <- apply(rbind(m1,m2),2,max)
# # de.gene <- names(max.m)[max.m > log2(11)]
# 
# 
# up.gene <- intersect(up.gene,de.gene)
# dn.gene <- intersect(dn.gene,de.gene)
# 
# 
# 
# g.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = 　dn.gene)
# 
# write.csv(x=g.df,quote = FALSE,file='~/Desktop/g.df.csv')
# 
# 
# # cor.with.SPARCL1 <- cor(PRI.log2.fpkm.matrix[c(dn.gene),] %>% t,method='spearman')[SPARCL1,]
# # 
# # #cor.with.SPARCL1 <- cor(PRI.log2.fpkm.matrix[c(NET.SI.DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene),] %>% t,method='spearman')[SPARCL1,]
# # 
# # cor.with.SPARCL1[c.dn.gene] %>% boxplot
# # 
# # 
# # hh <- names(cor.with.SPARCL1)[cor.with.SPARCL1 >= 0.5]
# 
# 
# g.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = 　hh)
# 
# 
# plot(PRI.log2.fpkm.matrix[c(PRRX1,dn.gene[1]),] %>% t)
# 
# y <- PRI.log2.fpkm.matrix[PRRX1,]
# 
# dn.gene <- setdiff(dn.gene,PRRX1)
# r.value <- foreach(g = dn.gene,.combine='c') %do% {
#     df <- data.frame(x=PRI.log2.fpkm.matrix[g,],y)  
#     fit.rs <- lm(df,formula = y ~ x)
#     mad(fit.rs$residuals)
# }
# 
# names(r.value) <- dn.gene



#####################
# median.tpm.matrix <- read.gct(file='client-side/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
# get.gene.id <- function(x) {
#   strsplit(x=x,split='\\.') %>% unlist %>% head(1)
# }
# rownames(median.tpm.matrix) <- sapply(rownames(median.tpm.matrix),get.gene.id)
# 
# 
# g <- intersect(rownames(median.tpm.matrix),BRCA.Basal.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene )
# pheatmap(log2(median.tpm.matrix[g,]+1))
# 
# 
# g <- intersect(rownames(median.tpm.matrix),BRCA.LumB.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene )
# pheatmap(log2(median.tpm.matrix[g,]+1))
# 
# 
# g <- intersect(rownames(median.tpm.matrix),COAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene )
# pheatmap(log2(median.tpm.matrix[g,]+1))
# 
# 
# g <- intersect(rownames(median.tpm.matrix),PRAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene )
# pheatmap(log2(median.tpm.matrix[g,]+1))




