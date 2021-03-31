# Aim: indepedent validation to show that M.vs.P DE analysis results could be confounded by sample biopsy site. 

# This script performs DE analysis between primary and liver-metastasis for breast cancer.Dataset is downloaded from SRA SRP043470.
# We show that liver-specific genes are identified as up-regulated in the comparision, suggeting it is necessary
# to remove them in the DE final gene list.

require(plyr)
require(dplyr)
require(DESeq2)
source('client-side/code/util.R')
load('server-side/RData/BRACA_SRP043470.RData')

p.sample               <- BRACA_SRP043470_Metadata$Run[BRACA_SRP043470_Metadata$source_name=='Primary Breast Tumour'] %>% as.character()
m.sample               <- BRACA_SRP043470_Metadata$Run[BRACA_SRP043470_Metadata$source_name=='Liver Metastasis Tumour'] %>% as.character()
protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character


g1             <- get.expressed.gene(BRACA_SRP043470_log2.fpkm.matrix[,p.sample])
g2             <- get.expressed.gene(BRACA_SRP043470_log2.fpkm.matrix[,m.sample])
expressed.gene <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)

expr.matrix    <- cbind(BRACA_SRP043470_log2.read.count.matrix[expressed.gene,p.sample],BRACA_SRP043470_log2.read.count.matrix[expressed.gene,m.sample])
expr.matrix    <- 2^expr.matrix - 1



df             <- data.frame(condition=c(rep(x='Primary',times=length(p.sample)), rep(x='Metastatic',times=length(m.sample))))
df$condition   <- factor(df$condition,levels = c('Primary','Metastatic'))

dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
                                         colData = df,
                                         design = ~ condition )
dds     <- DESeq(dds)
res     <- results(dds,contrast = c('condition','Metastatic','Primary')) %>% as.data.frame
res     <- res[order(res$pvalue),]
res     <- res[complete.cases(res),]     
up.gene <- rownames(res)[res$log2FoldChange > 1  & res$padj < 0.05]
dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]
SRP043470.de.res.metastasis.liver.vs.breast   <- res




# load('server-side//RData//Liver.RData')
# Ref.liver.log2.read.count.matrix <- log2.read.count.matrix
# Ref.liver.log2.fpkm.matrix       <- log2.fpkm.matrix
# liver.expressed.gene              <- get.expressed.gene(Ref.liver.log2.fpkm.matrix)


load('server-side//RData//Liver.RData')
female.sample                     <- sample.meta.df$sample.id[sample.meta.df$gender == 'Female'] %>% as.character()
Ref.liver.log2.read.count.matrix  <- log2.read.count.matrix[,female.sample]
Ref.liver.log2.fpkm.matrix        <- log2.fpkm.matrix[,female.sample]
liver.expressed.gene              <- get.expressed.gene(Ref.liver.log2.fpkm.matrix)

# load('server-side/RData/SRP068976_PairNor.RData')
# Ref.liver.log2.read.count.matrix  <- SRP068976_log2.read.count.matrix 
# Ref.liver.log2.fpkm.matrix        <- SRP068976_log2.fpkm.matrix 
# liver.expressed.gene              <- get.expressed.gene(Ref.liver.log2.fpkm.matrix)
# 


expressed.gene <- intersect(protein.coding.gene.id,c(liver.expressed.gene,g1) %>% unique)
expr.matrix    <- cbind(Ref.liver.log2.read.count.matrix[expressed.gene,],BRACA_SRP043470_log2.read.count.matrix[expressed.gene,p.sample])
expr.matrix    <- 2^expr.matrix - 1
df             <- data.frame(condition=c(rep(x='Liver',times=ncol(Ref.liver.log2.read.count.matrix)), rep(x='Primary',times=length(p.sample))))
df$condition   <- factor(df$condition,levels = c('Primary','Liver'))
dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
                                         colData = df,
                                         design = ~ condition )
dds     <- DESeq(dds)
res     <- results(dds,contrast = c('condition','Liver','Primary')) %>% as.data.frame
res     <- res[order(res$pvalue),]
res     <- res[complete.cases(res),]     
up.gene <- rownames(res)[res$log2FoldChange > 1  & res$padj < 0.05]
dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]
SRP043470.de.res.liver.vs.breast   <- res







save(file='client-side/output/validation.of.confounding.R.output/validation.of.confounding.RData',list=c('SRP043470.de.res.metastasis.liver.vs.breast','SRP043470.de.res.liver.vs.breast'))

################################## Trash code below this line ###################
# load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
# i.gene <- intersect(rownames(de.res.liver.vs.breast),rownames(SRP043470.liver.metastasis.vs.primary.res))
# plot(x=de.res.liver.vs.breast[i.gene,'log2FoldChange'],y=SRP043470.liver.metastasis.vs.primary.res[i.gene,'log2FoldChange'],xlim=c(-15,15),ylim=c(-15,15),xlab='Liver.vs.Breast',ylab='Metastsis.vs.Primary')
# lines(c(-15,15),c(-15,15))
# lowess(x=de.res.liver.vs.breast[i.gene,'log2FoldChange'],y=SRP043470.liver.metastasis.vs.primary.res[i.gene,'log2FoldChange']) %>% lines(lwd=5,col='red')
# 
# # Actually, here we could use Piecewise fitting to fit two lines:
# # y=0 & y= k.x + b.  For the genes that show high up-regulation, if in the M.vs.P analysis it up-regulates too much, then we still think it is a DE gene, instead of removing it!
# 
# 
# require(org.Hs.eg.db)
# up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)
