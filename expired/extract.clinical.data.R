library(cgdsr)
library(survival)
library(survminer)

mycgds          <-  CGDS("http://www.cbioportal.org/public-portal/")
study.id        <- 'brca_tcga'
mycaselist      <- getCaseLists(mycgds,study.id)
cast.list.id    <- 'brca_tcga_rna_seq_v2_mrna'
myclinicaldata  <- getClinicalData(mycgds,cast.list.id)

flag <- myclinicaldata$METASTATIC_SITE !='' 
data <- myclinicaldata[flag,]
id   <- rownames(data)
id   <- id[id != 'TCGA.BH.A18V.06']
conver.id <- function(x){
    tmp <- strsplit(x = x,split = "\\.") %>% unlist 
    paste(tmp,collapse  = '-')
  
}
id <- sapply(id,conver.id)

met.id <- id
pri.id <- setdiff(colnames(TCGA.breast.cancer.log2.fpkm.matrix),met.id)[1:24]

expr.matrix  <- cbind(TCGA.breast.cancer.log2.read.count.matrix[,met.id],TCGA.breast.cancer.log2.read.count.matrix[,pri.id])
expr.matrix  <- 2^expr.matrix - 1
flag         <- apply(expr.matrix,1,function(x) sum(x>=1)) # filter lowly expressed genes
expr.matrix  <- expr.matrix[flag == ncol(expr.matrix),]
df           <- data.frame(condition=c(rep(x='MET500',times=length(met.id)), rep(x='TCGA',times=length(pri.id))))
df$condition <- factor(df$condition,levels = c('TCGA','MET500'))
dds          <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
                                       colData = df,
                                       design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$pvalue),]
res <- res[complete.cases(res),]



up.de.gene   <- rownames(res)[res$padj < 0.001 & res$log2FoldChange >  1.5]
down.de.gene <- rownames(res)[res$padj < 0.001 & res$log2FoldChange < -1.5]





gene <- c(liver.metastasis.vs.primary.tumor.non.Basal.up.de.gene,liver.metastasis.vs.primary.tumor.non.Basal.down.de.gene)
col.vec <- ifelse(colnames(TCGA.breast.cancer.log2.fpkm.matrix) %in% id,'red','black')
pca.rs  <- prcomp(TCGA.breast.cancer.log2.fpkm.matrix[gene,] %>% t)
plot(pca.rs$x[,1:2],pch=19,col=col.vec)
