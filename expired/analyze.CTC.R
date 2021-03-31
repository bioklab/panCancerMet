load('server-side/RData/GSE111842.RData')
#log2.fpkm.matrix <- log2.fpkm.matrix[,colnames(log2.fpkm.matrix) %in% c('GSM3040933','GSM3040934','GSM3040935') == FALSE]
#GSE111842_Metadata <- read.csv("~/Project/BreastCancerMetaPotenial/server-side/RData/GSE111842_Metadata.csv", stringsAsFactors=FALSE)
#GSE111842_Metadata <- GSE111842_Metadata[rownames(GSE111842_Metadata) %in% c('GSM3040933','GSM3040934','GSM3040935') == FALSE,]
GSE111842_Metadata <- read.csv("~/Project/BreastCancerMetaPotenial/server-side/RData/GSE111842_Metadata.csv", stringsAsFactors=FALSE)
rownames(GSE111842_Metadata) <- GSE111842_Metadata$Sample_Name


non.zero.gene.cnt <- apply(log2.fpkm.matrix,2,function(x) sum(x>0))
non.zero.gene.cnt <- sort(non.zero.gene.cnt,decreasing = TRUE)
cbind(GSE111842_Metadata[names(non.zero.gene.cnt),],non.zero.gene.cnt)


ctc.sample           <- GSE111842_Metadata$Sample_Name[GSE111842_Metadata$sample_type == 'circulating tumor cells sorted from peripheral blood']
ctc.log2.fpkm.matrix <- log2.fpkm.matrix[,ctc.sample]
ctc.log2.tpm.matrix  <- log2.tpm.matrix[,ctc.sample]




dist.obj       <- as.dist(1- cor(dd[pam50.gene,],method='spearman'))

tsne.rs        <- Rtsne(dist.obj,perplexity = 5)
#tsne.rs        <- Rtsne(log2.fpkm.matrix %>% t,perplexity = 8,check_duplicates = FALSE)

rownames(tsne.rs$Y) <- colnames(dd)
plot(tsne.rs$Y)



pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') # well, I borrow some information from the metastatic breast cancer evaluation project 
pam50.gene                  <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- ctc.log2.fpkm.matrix[pam50.gene.df$probe %>% as.character,] %>% t 
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )




pca.rs <- prcomp(log2.fpkm.matrix %>% t)
plot(pca.rs$x[,1:2])



pc1 <- pca.rs$x[,1]
outlier.sample <- names(pc1)[pc1 < -200]

log2.fpkm.matrix <- log2.fpkm.matrix[,colnames(log2.fpkm.matrix) %in% outlier.sample == FALSE]
GSE111842_Metadata <- GSE111842_Metadata[rownames(GSE111842_Metadata) %in% outlier.sample == FALSE,]

flag <- apply(log2.fpkm.matrix,1,function(x) median(x) > log2(1+0.1))

pca.rs <- prcomp(log2.fpkm.matrix[flag,] %>% t)
plot(pca.rs$x[,1:2])


draw.df <- data.frame(pca.rs$x[rownames(GSE111842_Metadata),1:10],cell.type=GSE111842_Metadata$sample_type)
ggplot(draw.df) + geom_point(aes(x=PC1,y=PC2,color=cell.type))


draw.df$TMSB4X <- log2.fpkm.matrix['ENSG00000205542',rownames(draw.df)]
draw.df$DONSON <- log2.fpkm.matrix['ENSG00000159147',rownames(draw.df)]
draw.df$TRIM9 <- log2.fpkm.matrix['ENSG00000100505',rownames(draw.df)]
draw.df$SOCS7 <- log2.fpkm.matrix['ENSG00000205176',rownames(draw.df)]

ggplot(draw.df) + geom_boxplot(aes(y=SOCS7,color=cell.type))





r1 <- pca.rs$rotation[,'PC1'] %>% sort
head(r1,10)
tail(r1,10)




load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
g <- intersect(CCLE.rna.seq.marker.gene.1000,rownames(log2.fpkm.matrix))


dist.obj       <- as.dist(1- cor(log2.fpkm.matrix[g,],method='spearman'))
dist.obj       <- as.dist(1- cor(log2.fpkm.matrix[pam50.gene,],method='spearman'))

set.seed(8) # I want to reproduce the tsne results, 8 is just a arbitrary numnber, it DOES NOT change the conclusion
tsne.rs        <- Rtsne(dist.obj,perplexity = 8)
#tsne.rs        <- Rtsne(log2.fpkm.matrix %>% t,perplexity = 8,check_duplicates = FALSE)

rownames(tsne.rs$Y) <- colnames(log2.fpkm.matrix)

draw.df <- data.frame(pca.rs$x[rownames(GSE111842_Metadata),1:5],tsne.rs$Y[rownames(GSE111842_Metadata),],cell.type=GSE111842_Metadata$sample_type)
ggplot(draw.df) + geom_point(aes(x=X1,y=X2,color=cell.type))
ggplot(draw.df) + geom_point(aes(x=PC2,y=PC3,color=cell.type))
draw.df$TMSB4X <- log2.fpkm.matrix['ENSG00000205542',rownames(draw.df)]
draw.df$DONSON <- log2.fpkm.matrix['ENSG00000159147',rownames(draw.df)]
draw.df$TRIM9 <- log2.fpkm.matrix['ENSG00000100505',rownames(draw.df)]

ggplot(draw.df) + geom_boxplot(aes(y=TRIM9,color=cell.type))

