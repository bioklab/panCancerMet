source('client-side/code/DEBoost.R')

################################################################################################################
# Prepare liver data, here we only want male samples
################################################################################################################
load('server-side//RData//Liver.RData')
Male.sample                   <- sample.meta.df$sample.id[sample.meta.df$gender == 'Male'] %>% as.character()
REF.log2.read.count.matrix    <- log2.read.count.matrix[,Male.sample]
REF.log2.tpm.matrix           <- log2.tpm.matrix[,Male.sample]


load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')

################################################################################################################
# Prepare primary cancer data
################################################################################################################
load('server-side/RData//Prostate Adenocarcinoma.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.prostate.cancer.sample]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]


################################################################################################################
# Prepare metastatic cancer data
################################################################################################################
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
MET.log2.read.count.matrix <- MET500.log2.read.count.matrix[,pure.MET.prostate.cancer.sample]
MET.log2.tpm.matrix       <- MET500.log2.tpm.matrix[,pure.MET.prostate.cancer.sample]



PRAD.DE.rs <- perform.DE.analysis.between.primary.and.metastatic.cancer(
  PRI.log2.tpm.matrix = PRI.log2.tpm.matrix, PRI.log2.read.count.matrix = PRI.log2.read.count.matrix,
  MET.log2.tpm.matrix = MET.log2.tpm.matrix, MET.log2.read.count.matrix = MET.log2.read.count.matrix,
  REF.log2.tpm.matrix = REF.log2.tpm.matrix, REF.log2.read.count.matrix = REF.log2.read.count.matrix,
  TCGA.best.cell.line = TCGA.best.cell.line, MET500.best.cell.line = MET500.best.cell.line
)


save(file='client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData',list=c('PRAD.DE.rs'))








############# Trash code #####################

# up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = DE.rs$tumor.intrinsic.DE.gene.rs$up.gene) 
# GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
# GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 & GO.rs.1.up$pvalue < 0.01, ]
# 
# dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene) 
# GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
# GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 & GO.rs.1.dn$pvalue < 0.01, ]
# 
# 

#######################################################################################################################################################################################
###### Prepare the data 
#######################################################################################################################################################################################

# MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# 
# #######################################################################################################################################################################################
# ###### Function to perform DE analysis between tumor samples
# #######################################################################################################################################################################################
# perform.DE.analysis.between.metastatic.and.primary.cancer <- function(){
#   g1             <- get.expressed.gene(TCGA.prostate.cancer.log2.tpm.matrix[,TCGA.sample])
#   g2             <- get.expressed.gene(MET500.log2.tpm.matrix[,MET500.sample])
#   expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
#   
#   expr.matrix    <- cbind(MET500.log2.read.count.matrix[expressed.gene,MET500.sample],TCGA.prostate.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample])
#   expr.matrix    <- 2^expr.matrix - 1
#   
#   
#   
#   df             <- data.frame(condition=c(rep(x='MET500',times=length(MET500.sample)), rep(x='TCGA',times=length(TCGA.sample))))
#   df$condition   <- factor(df$condition,levels = c('TCGA','MET500'))
#   
#   dds           <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                           colData = df,
#                                           design = ~ condition)
#   dds <- DESeq(dds)
#   res <- results(dds,contrast = c('condition','MET500','TCGA')) %>% as.data.frame
#   res <- res[order(res$pvalue),]
#   res <- res[complete.cases(res),]     
#   res
# }
# 
# perform.DE.analysis.between.ref.tissue.and.primary.cancer <- function(){
#   g1             <- get.expressed.gene(TCGA.prostate.cancer.log2.tpm.matrix[,TCGA.sample])
#   g2             <- get.expressed.gene(Ref.liver.log2.tpm.matrix)
#   expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
#   
#   expr.matrix    <- cbind(Ref.liver.log2.read.count.matrix[expressed.gene,],TCGA.prostate.cancer.log2.read.count.matrix[expressed.gene,TCGA.sample])
#   expr.matrix    <- 2^expr.matrix - 1
#   
#   df             <- data.frame(condition=c(rep(x='ref.tissue',times=ncol(Ref.liver.log2.read.count.matrix)), rep(x='TCGA',times=length(TCGA.sample))))
#   df$condition   <- factor(df$condition,levels = c('TCGA','ref.tissue'))
#   
#   dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                            colData = df,
#                                            design = ~ condition )
#   dds <- DESeq(dds)
#   res <- results(dds,contrast = c('condition','ref.tissue','TCGA')) %>% as.data.frame
#   res <- res[order(res$pvalue),]
#   res <- res[complete.cases(res),]     
#   res
# }
# 
# 
# 
# 
# 
# DEBoost.filtering <- function(){
#   up.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange > 1  & deseq2.res$padj < 0.05]
#   dn.gene            <- rownames(deseq2.res)[deseq2.res$log2FoldChange < -1 & deseq2.res$padj < 0.05]
#   
#   ######### piecewise linear regression to identify ref-tissue specific genes #########
#   c.gene        <- intersect(rownames(deseq2.ref.res),rownames(deseq2.res))
#   x             <- deseq2.ref.res[c.gene,'log2FoldChange']
#   y             <- deseq2.res[c.gene,'log2FoldChange']
#   lin.mod       <- lm(y~x)
#   segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)    
#   tmp           <- summary(segmented.mod)
#   psi           <- tmp$psi[1,'Est.']
#   ref.specific.gene  <- rownames(deseq2.ref.res)[deseq2.ref.res$log2FoldChange > psi]
#   
#   ######### gene filtering based on correlation with purity#########
#   # p              <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
#   # cor.vec        <- cor(TCGA.breast.cancer.log2.tpm.matrix[,TCGA.sample] %>% t,p,method='spearman') %>% c
#   # names(cor.vec) <- rownames(TCGA.breast.cancer.log2.tpm.matrix)
#   # cor.data       <- cor.vec[c(up.gene,dn.gene)]
#   # model          <- normalmixEM(cor.data)
#   # mu.vec         <- model$mu
#   # sigma.vec      <- model$sigma
#   # if(mu.vec[1] < mu.vec[2]){
#   #     mu     <- mu.vec[1]
#   #     sigma  <- sigma.vec[1]
#   # } else {
#   #     mu      <- mu.vec[2]
#   #     sigma   <- sigma.vec[2]
#   # }
#   # cor.cut.off  <- mu + 3.0 * sigma
#   # neg.cor.gene <- names(cor.data)[ cor.data <= cor.cut.off]
#   # ggplot.graph <- ggplot(data.frame(x=cor.data), aes(x=x)) + geom_histogram() + geom_density(aes(y=..density.. * 50)) + geom_vline(xintercept = cor.cut.off)
#   # 
#   
#   ####### gene filtering based on expression level #############
#   immune.gene.list        <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/immune.gene.list.csv", stringsAsFactors=FALSE)$x
#   TCGA.median             <- apply(TCGA.prostate.cancer.log2.tpm.matrix[, TCGA.sample], 1, median)
#   MET500.median           <- apply(MET500.log2.tpm.matrix[, MET500.sample], 1, median)
#   median.max              <- apply(rbind(TCGA.median,MET500.median),2,max)
#   x                       <- median.max[immune.gene.list]
#   x                       <- x[is.na(x) == FALSE]
#   q                       <- quantile(x)
#   tpm.cut.off            <- q['75%'] + 1.5 * IQR(x)
#   low.expr.gene           <- names(median.max)[median.max < tpm.cut.off]
#   
#   
#   ####### gene filtering based on wilcoxon rank test #######################
#   wilcox.test.p.value.vec <- foreach(g=c(up.gene,dn.gene),.combine='c') %do% {
#     wilcox.test(MET500.log2.tpm.matrix[g,MET500.sample],TCGA.prostate.cancer.log2.tpm.matrix[g,TCGA.sample])$p.value  
#   }
#   names(wilcox.test.p.value.vec) <- c(up.gene,dn.gene)
#   not.robust.gene                <- names(wilcox.test.p.value.vec) [ wilcox.test.p.value.vec > 0.05]
#   
#   
#   uDE.gene <- c(ref.specific.gene,low.expr.gene,not.robust.gene) %>% unique()
#   up.gene  <- setdiff(up.gene,uDE.gene)
#   dn.gene  <- setdiff(dn.gene,uDE.gene)
#   
#   
#   list(up.gene=up.gene,dn.gene=dn.gene)
# }
# 
# 
# 
# #######################################################################################################################################################################################
# ###### DE analysis
# #######################################################################################################################################################################################
# 
# TCGA.sample    <- pure.TCGA.prostate.cancer.polyA.sample
# 
# 
# 
# deseq2.res      <- perform.DE.analysis.between.metastatic.and.primary.cancer()
# deseq2.ref.res  <- perform.DE.analysis.between.ref.tissue.and.primary.cancer()
# de.gene.rs      <- DEBoost.filtering()
# # up.gene         <- de.gene.rs$up.gene
# # dn.gene         <- de.gene.rs$dn.gene
# 
# # write.csv(x=basal.up.gene,file = 'client-side/output/DE.breast.cancer.R.output/basal.up.csv',quote=FALSE)
# # write.csv(x=basal.dn.gene,file = 'client-side/output/DE.breast.cancer.R.output/basal.dn.csv',quote=FALSE)
# 
# up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = de.gene.rs$up.gene) 
# GO.rs.1.up             <- enrichGO(gene=up.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
# GO.rs.1.up             <- GO.rs.1.up[ GO.rs.1.up$Count >= 5 & GO.rs.1.up$pvalue < 0.01, ]
# 
# 
# dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = de.gene.rs$dn.gene) 
# GO.rs.1.dn             <- enrichGO(gene=dn.gene.annotation.df$SYMBOL,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
# GO.rs.1.dn             <- GO.rs.1.dn[ GO.rs.1.dn$Count >= 5 & GO.rs.1.dn$pvalue < 0.01, ]
# 
# 
# 
# G1.and.S <- read.table("~/Project/BreastCancerMetaPotenial/client-side/Data/cell.cycle/G1.and.S.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
# G2.and.M <- read.table("~/Project/BreastCancerMetaPotenial/client-side/Data/cell.cycle/G2.and.M.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
# (intersect(G1.and.S,up.gene.annotation.df$SYMBOL) %>% length) / (G1.and.S %>% length)
# (intersect(G2.and.M,up.gene.annotation.df$SYMBOL) %>% length) / (G2.and.M %>% length)
# 
# 
# G1.and.S.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'SYMBOL',columns = c('ENSEMBL','SYMBOL'),keys = G1.and.S) 
# G1.and.S.es  <- intersect(G1.and.S.df$ENSEMBL,rownames(TCGA.prostate.cancer.log2.tpm.matrix))
# 
# 
# G2.and.M.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'SYMBOL',columns = c('ENSEMBL','SYMBOL'),keys = G2.and.M) 
# G2.and.M.es  <- intersect(G2.and.M.df$ENSEMBL,rownames(TCGA.prostate.cancer.log2.tpm.matrix))
# 
# 
# tcga <- apply(TCGA.prostate.cancer.log2.tpm.matrix[G1.and.S.es,],1,median)
# met  <- apply(MET500.log2.tpm.matrix[G1.and.S.es,MET500.sample],1,median)
# boxplot(tcga,met,main='G1/S')
# wilcox.test(tcga,met)
# 
# tcga <- apply(TCGA.prostate.cancer.log2.tpm.matrix[G2.and.M.es,],1,median)
# met  <- apply(MET500.log2.tpm.matrix[G2.and.M.es,MET500.sample],1,median)
# boxplot(tcga,met,main='G2/M')
# wilcox.test(tcga,met)
# 
# 
# 
# 
# 
# 
# 
# 
# # get.spectral.clustering.coordinates <- function(A,no.of.eigen.vector=0){ 
# #   diag(A)  <-   0
# #   A        <-  abs(A)
# #   d        <-  apply(A,1,sum)
# #   I        <-  diag(1,nrow = nrow(A)) 
# #   L        <-   I -  diag(1/sqrt(d)) %*% A %*% diag(1/sqrt(d))
# #   
# #   
# #   #L is positive semi-definite and have n non-negative real-valued eigenvalues
# #   tmp             <-  eigen(L)
# #   eigen.values    <-  tmp$values[length(tmp$values):1]
# #   eigen.vectors   <-  tmp$vectors[,length(tmp$values):1]
# #   
# #   #automaticaly choose number of eigen-vectors
# #   q1              <-  quantile(eigen.values)[2]
# #   q3              <-  quantile(eigen.values)[4]
# #   low             <-  q1 - 3*(q3 - q1)
# #   up              <-  q1 + 3*(q3 - q1)
# #   k               <-  sum(eigen.values < low) # pick out eigen-values that are relatively small
# #   
# #   if(no.of.eigen.vector > 1){
# #     k <- no.of.eigen.vector
# #   }
# #   if(k == 1){
# #     list(eigen.values=eigen.values,no.of.eigen.vector=1)
# #   }else{
# #     coordinate.matrix            <- eigen.vectors[,2:(k+1)] # pick out the eigen-vectors associated with the k smallest non-zero eigen-values
# #     rownames(coordinate.matrix)  <- rownames(A)
# #     colnames(coordinate.matrix)  <- paste('eigen',1:(k),sep="-")
# #     list(coordinates=coordinate.matrix,eigen.values=eigen.values,no.of.eigen.vector=k)
# #   }
# # }
# # 
# # 
# # dn.gene.co.expr.matrix <- cor(TCGA.prostate.cancer.log2.tpm.matrix[de.gene.rs$up.gene,] %>% t,method='spearman')
# # 
# # 
# # rs <- get.spectral.clustering.coordinates(dn.gene.co.expr.matrix %>% abs)
# # 
# # x <- rs$coordinates[,1] %>% sort
# # pheatmap(abs(dn.gene.co.expr.matrix[names(x),names(x)]),cluster_rows=F,cluster_cols=F)
