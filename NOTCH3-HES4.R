#Clonality of circulating tumor cells in breast cancer brain metastasis patients. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6720990/#CR21
#Notch3 Maintains Luminal Phenotype and Suppresses Tumorigenesis and Metastasis of Breast Cancer via Trans-Activating Estrogen Receptor-α. NOTCH3, repress EMT and metastasis
#Constitutive NOTCH3 signaling promotes the growth of basal breast cancers. https://cancerres.aacrjournals.org/content/77/6/1439.long
#Delta-notch--and then? Protein interactions and proposed modes of repression by hes and hey bhlh factors
################################# cancer testis antigen , APOE, MYC,breast cancer pathway gene (cutoff 0.5) , MMP23B, FZD8 wnt pathway, TNFRSF18(GITR, Treg) immune checkpoint ##########
# Ref: Rationale for anti-GITR cancer immunotherapy
# Ref:SPIB and BATF provide alternate determinants of IRF4 occupancy in diffuse large B-cell lymphoma linked to disease heterogeneity 
# The Hes gene family: repressors and oscillators that orchestrate embryogenesis. HES4 exists in human not in mouse
###########################################
source('client-side/code/util.R')
require(DESeq2)
require(org.Hs.eg.db)
require(plyr)
require(dplyr)


get.expressed.gene <- function(expr.matrix,cut.off=1){
  m.expr <- apply(expr.matrix,1,median)
  rownames(expr.matrix)[m.expr >= cut.off ]
}

############ DE between HES4-high vs HES4-low samples, TCGA samples ###########
load('server-side/RData/Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')

TCGA.breast.cancer.log2.tpm.matrix        <- log2.tpm.matrix
TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix

HES4       <- 'ENSG00000188290'
NOTCH1     <- 'ENSG00000148400'
NOTCH3     <- 'ENSG00000074181'
HES4.expr  <- TCGA.breast.cancer.log2.tpm.matrix[HES4,pure.PRI.breast.cancer.Basal.sample] %>% sort
hist(HES4.expr,breaks=30)
quantile(HES4.expr,probs = seq(0, 1, 0.05))  %>% plot
# plot(TCGA.breast.cancer.log2.tpm.matrix[c(NOTCH3,HES4),] %>% t)
# cor(TCGA.breast.cancer.log2.tpm.matrix[c(NOTCH3,HES4),] %>% t,method='spearman')
# 
# abline(h=3)


l.sample      <- names(HES4.expr)[HES4.expr <  median(HES4.expr)]
h.sample      <- names(HES4.expr)[HES4.expr >=  median(HES4.expr)]
l.expr.matrix <- TCGA.breast.cancer.log2.tpm.matrix[,l.sample]
h.expr.matrix <- TCGA.breast.cancer.log2.tpm.matrix[,h.sample]

g1                     <- get.expressed.gene(l.expr.matrix)
g2                     <- get.expressed.gene(h.expr.matrix)
expressed.gene         <- c(g1,g2) %>% unique

tmp                  <- cbind(TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,h.sample],TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,l.sample])
read.count.matrix    <- 2^tmp - 1


df             <- data.frame(condition=c(rep(x='HES4.high',times=length(h.sample)), rep(x='HES4.low',times=length(l.sample))))
df$condition   <- factor(df$condition,levels = c('HES4.low','HES4.high'))

dds            <- DESeqDataSetFromMatrix(countData = round(read.count.matrix),
                                         colData = df,
                                         design = ~ condition )
dds            <- DESeq(dds)
res            <- results(dds,contrast = c('condition','HES4.high','HES4.low')) %>% as.data.frame
res            <- res[order(res$pvalue),]
res            <- res[complete.cases(res),]     
res            <- as.data.frame(res)


HES4.h.up.gene    <- rownames(res)[res$log2FoldChange >  1 & res$padj < 0.05]
HES4.h.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = HES4.h.up.gene)
HES4.h.dn.gene    <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]
HES4.h.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = HES4.h.dn.gene)

TCGA.HES4.high.vs.low.res <- res 



############ DE between HES4-high vs HES4-low samples, SRP157974 samples ###########
load('server-side/RData/SRP157974_PrimaryTumor.RData')
require(Rtsne)
require(genefu)
TCGA.breast.cancer.polyA.sample <- sample.meta.df$sample.id[sample.meta.df$primary.disease.or.tissue == 'Breast Invasive Carcinoma']

#### Perform pam50 subtyping for TCGA polyA samples###########
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') # well, I borrow some information from the metastatic breast cancer evaluation project 
pam50.gene                  <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- cbind(TCGA.breast.cancer.log2.tpm.matrix[pam50.gene.df$probe %>% as.character,TCGA.breast.cancer.polyA.sample],
                                     SRP157974_PrimaryTumor_log2.tpm.matrix[pam50.gene.df$probe %>% as.character,]
) %>% t
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )


dist.obj       <- as.dist(1- cor(pam50.gene.expr %>% t,method='spearman'))
set.seed(8) # I want to reproduce the tsne results, 8 is just a arbitrary numnber, it DOES NOT change the conclusion
tsne.rs        <- Rtsne(dist.obj,perplexity = 15)

draw.df <- data.frame(dim1=tsne.rs$Y[,1],
                      dim2=tsne.rs$Y[,2],
                      subtype       = pam50.subtype.rs$subtype[rownames(pam50.gene.expr)]
)
ggplot(draw.df,aes(x=dim1,y=dim2,color=subtype)) + geom_point(size=6) +  scale_shape_manual(values=c(8,15:18))

basal.sample <- names(pam50.subtype.rs$subtype)[pam50.subtype.rs$subtype == 'Basal']
basal.sample <- basal.sample[grepl(x=basal.sample,pattern = 'TCGA') == FALSE]



load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000   <- names(sort(rank.sd,decreasing =TRUE))[1:1000]

rs <- pick.out.cell.line(expr.of.samples = SRP157974_PrimaryTumor_log2.fpkm.matrix[,basal.sample], expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
r  <- apply(rs$correlation.matrix,1,function(x) rank(x)['HDQP1_BREAST'])
plot(sort(r))
basal.sample <- names(r)[r >= 1000]


cDNA.sample  <- SRP157974_PrimaryTumor_Metadata$Run[SRP157974_PrimaryTumor_Metadata$LibrarySelection =='cDNA'] %>% as.character()
basal.sample <- intersect(basal.sample,cDNA.sample)

SRP157974_PrimaryTumor_log2.tpm.matrix       <- SRP157974_PrimaryTumor_log2.tpm.matrix[,basal.sample]
SRP157974_PrimaryTumor_log2.read.count.matrix <- SRP157974_PrimaryTumor_log2.read.count.matrix[,basal.sample]


HES4.expr <- SRP157974_PrimaryTumor_log2.tpm.matrix[HES4,] %>% sort
hist(HES4.expr,breaks=40)
quantile(HES4.expr,probs = seq(0, 1, 0.05))  %>% plot


l.sample      <- names(HES4.expr)[HES4.expr  <  median(HES4.expr)]
h.sample      <- names(HES4.expr)[HES4.expr >=  median(HES4.expr)]
l.expr.matrix <- SRP157974_PrimaryTumor_log2.tpm.matrix[,l.sample]
h.expr.matrix <- SRP157974_PrimaryTumor_log2.tpm.matrix[,h.sample]

g1                     <- get.expressed.gene(l.expr.matrix)
g2                     <- get.expressed.gene(h.expr.matrix)
expressed.gene         <- c(g1,g2) %>% unique

tmp    <- cbind(SRP157974_PrimaryTumor_log2.read.count.matrix[expressed.gene,h.sample],SRP157974_PrimaryTumor_log2.read.count.matrix[expressed.gene,l.sample])
read.count.matrix    <- 2^tmp - 1
df             <- data.frame(condition=c(rep(x='HES4.high',times=length(h.sample)), rep(x='HES4.low',times=length(l.sample))))
df$condition   <- factor(df$condition,levels = c('HES4.low','HES4.high'))
dds            <- DESeqDataSetFromMatrix(countData = round(read.count.matrix),
                                         colData = df,
                                         design = ~ condition  )
dds            <- DESeq(dds)
res            <- results(dds,contrast = c('condition','HES4.high','HES4.low')) %>% as.data.frame
res            <- res[order(res$pvalue),]
res            <- res[complete.cases(res),]     
res            <- as.data.frame(res)
SRP157974.HES4.high.vs.low.res <- res 
SRP157974.HES4.high.vs.low.res[NOTCH3,]

save(file='client-side/output/NOTCH3-HES4.R.output/NOTCH3-HES4.RData',list =c ('SRP157974.HES4.high.vs.low.res','TCGA.HES4.high.vs.low.res'))
















############ Compare NOTCH3  CNV between primary and met, basal-like subtype ###########
# load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
# load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
# load('client-side/output/TCGA.breast.cancer.meta.R.output/TCGA.breast.cancer.meta.RData')
# load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')
# load('client-side/output/organize.TCGA.and.MET500.breast.cancer.cnv.data.R.output/organize.TCGA.and.MET500.breast.cancer.cnv.data.RData')
# 
# MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# MET500.sample        <- MET500.breast.cancer.polyA.Basal.sample
# MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Basal.sample
# 
# MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
# TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))
# 
# boxplot(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
# wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
# 
# #boxplot(MET500.breast.cancer.cnv.matrix['NOTCH1',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH1',TCGA.subject.id])
# #wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH1',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH1',TCGA.subject.id])
# 
# 
# 
# ########## Compare CNV profile between primary and metastatic breast cancer
# gene.vec <- intersect(rownames(MET500.breast.cancer.cnv.matrix), rownames(TCGA.breast.cancer.cnv.matrix))
# df <- foreach(g = gene.vec,.combine='rbind') %do% {
#   df1                  <- data.frame(cnv= MET500.breast.cancer.cnv.matrix[g,MET500.subject.id], site= rep('MET500',length(MET500.subject.id)))
#   df2                  <- data.frame(cnv= TCGA.breast.cancer.cnv.matrix  [g,TCGA.subject.id],   site= rep('TCGA',length(TCGA.subject.id)))
#   draw.df              <- rbind(df1,df2)
#   draw.df$site         <- factor(draw.df$site,levels = c('TCGA','MET500'))
#   
#   #ggplot(draw.df) + geom_boxplot(aes(x=site,y=cnv),outlier.shape=NA,lwd=3) + ggplot.style + geom_jitter(aes(x=site,y=cnv),size=5.5) + xlab('')
#   p.value <- wilcox.test(MET500.breast.cancer.cnv.matrix[g,MET500.subject.id],TCGA.breast.cancer.cnv.matrix[g,TCGA.subject.id])$p.value
#   diff    <- median(MET500.breast.cancer.cnv.matrix[g,MET500.subject.id]) - median(TCGA.breast.cancer.cnv.matrix[g,TCGA.subject.id])
#   data.frame(p.value,diff,met=median(MET500.breast.cancer.cnv.matrix[g,MET500.subject.id]))
#   
#   #ggplot(draw.df) + geom_boxplot(aes(x=site,y=cnv),outlier.shape=NA,lwd=3) + ggplot.style + geom_jitter(aes(x=site,y=cnv),size=5.5) + xlab('')
#   
# }
# rownames(df) <- gene.vec
# df <- df[complete.cases(df),]
# df$fdr <- p.adjust(df$p.value,method='fdr')
# 
# 
# 
# 
# ############ Up-regulated genes are enirched at chr19 p13.12 ###########
# 
# chr19.p13.12.start <- 14000001
# chr19.p13.12.end   <- 16300000
# 
# flag              <- hg19.gene.info$chrom == '19' & hg19.gene.info$start > chr19.p13.12.start & hg19.gene.info$end < chr19.p13.12.end
# chr19.p13.12.gene <- hg19.gene.info[flag,'genename'] %>% as.character()
# flag              <- grepl(chr19.p13.12.gene,pattern = 'MIR') | grepl(chr19.p13.12.gene,pattern = 'LINC') | grepl(chr19.p13.12.gene,pattern = 'LOC') |grepl(chr19.p13.12.gene,pattern = 'SNORA')
# chr19.p13.12.gene <- chr19.p13.12.gene[!flag]
# 
# basal.up.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
# basal.dn.gene <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
# 
# c.gene           <- intersect(basal.up.gene,chr19.p13.12.gene)
# c.gene.df       <- hg19.gene.info[hg19.gene.info$genename %in% c.gene,]
# c.gene.df       <- c.gene.df[order(c.gene.df$start),]
# c.gene.df$chrom <- paste('chr',c.gene.df$chrom,sep='')
# write.table(x=c.gene.df[,c('chrom','start','end','genename')],file='client-side/output/NOTCH3-HES4.R.output/basal.up.chr19.p13.12.gene.bed',quote=FALSE,row.names=FALSE)
# 
# 
# # extract CNV segment call covering NOTHC3
# get.id <- function(x){
#     l <- strsplit(x = x,split='\\.')  %>% unlist 
#     l[1]
# }
# MET500.seg.data                         <- read.csv("~/Project/Cancer2CellLine/client-side/Data/CNV/cnv_v4.csv", stringsAsFactors=FALSE)
# MET500.seg.data                         <- MET500.seg.data[,c('Pipeline_ID','Chr','Start','End','Log2_Coverage_Ratio')]
# colnames(MET500.seg.data)               <- c('ID','chrom','loc.start','loc.end','seg.mean')
# MET500.seg.data$ID                      <- sapply(MET500.seg.data$ID,get.id)
# 
# 
# extract.segment.covering.gene <- function(gene.name) {
#     row     <-  hg19.gene.info[hg19.gene.info$genename == gene.name,]
#     chrom   <- row$chrom
#     start   <- row$start
#     end     <- row$end
#     flag    <- (MET500.seg.data$loc.start < start) & (MET500.seg.data$loc.end > end) & MET500.seg.data$ID %in% MET500.subject.id & MET500.seg.data$chrom == chrom
#     gene.seg        <- MET500.seg.data[flag,]
#     gene.seg$length <- gene.seg$loc.end - gene.seg$loc.start
#     gene.seg$chrom  <- paste('chr',gene.seg$chrom,sep = '')
#     gene.seg
# }
# NOTCH3.seg <- extract.segment.covering.gene('NOTCH3')
# write.table(NOTCH3.seg[,c('chrom','loc.start','loc.end')],sep = '\t',file='client-side/output/NOTCH3-HES4.R.output/NOTCH3.seg.bed',quote = FALSE,row.names=FALSE)
# 
# CHAC1.seg <- extract.segment.covering.gene('CHAC1')
# MAML1.seg <- extract.segment.covering.gene('MAML1')
# NRARP.seg <- extract.segment.covering.gene('NRARP')
# NOTCH1.seg <- extract.segment.covering.gene('NOTCH1')

# flag                                    <- MET500.seg.data$ID %in% MET500.subject.id & MET500.seg.data$chrom == '19'
# chr19.data                              <- MET500.seg.data[flag,]
# NOTCH3.start      <- hg19.gene.info[hg19.gene.info$genename == 'NOTCH3','start']
# NOTCH3.end        <- hg19.gene.info[hg19.gene.info$genename == 'NOTCH3','end']
# flag              <- (chr19.data$loc.start < NOTCH3.start) & (chr19.data$loc.end > NOTCH3.end)
# NOTCH3.seg        <- chr19.data[flag,]
# NOTCH3.seg$length <- NOTCH3.seg$loc.end - NOTCH3.seg$loc.start
# NOTCH3.seg$chrom  <- paste('chr',NOTCH3.seg$chrom,sep = '')









# ############# Compare NOTCH3 cnv between primary and met, LumB subtype ###########
# MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# MET500.sample        <- MET500.breast.cancer.polyA.LumB.sample
# MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Basal.sample
# 
# MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
# TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))
# 
# boxplot(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
# wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
# 
# 
# ############# Compare NOTCH3 cnv between primary and met, Her2 subtype ###########
# MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# MET500.sample        <- MET500.breast.cancer.polyA.Her2.sample
# MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Basal.sample
# 
# MET500.subject.id    <- MET500.sample.meta[MET500.sample %>% as.character(),'MET500.id']
# TCGA.subject.id      <- intersect(TCGA.sample,colnames(TCGA.breast.cancer.cnv.matrix))
# 
# boxplot(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])
# wilcox.test(MET500.breast.cancer.cnv.matrix['NOTCH3',MET500.subject.id],TCGA.breast.cancer.cnv.matrix['NOTCH3',TCGA.subject.id])


###############################################
save(file='client-side/output/NOTCH3-HES4.R.output/NOTCH3-HES4.RData',list=c('SRP157974.HES4.high.vs.low.res','TCGA.HES4.high.vs.low.res'))







############## analyze scRNA TNBC data ###########
# GSE118389_tpm_rsem <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/GSE118389_tpm_rsem.txt", row.names=1, stringsAsFactors=FALSE)
# GSE118389_tpm_rsem.matrix <- as.matrix(GSE118389_tpm_rsem)
# expr.vec <- GSE118389_tpm_rsem.matrix['HES4',]
# h.name <- names(expr.vec)[expr.vec > 0]
# l.sample <- names(expr.vec)[expr.vec  == 0]
# 
# p.value.vec <- foreach(g=rownames(GSE118389_tpm_rsem.matrix),.combine='c') %do% {
#     wilcox.test(GSE118389_tpm_rsem.matrix[g,h.name],GSE118389_tpm_rsem.matrix[g,l.sample])$p.value  
# }
# names(p.value.vec) <- rownames(GSE118389_tpm_rsem.matrix)
# 
# 
# log2.matrix <- log2(GSE118389_tpm_rsem.matrix + 1)
# 
# quantile(log2.matrix['NOTCH3',l.sample],probs=seq(0,1,0.05))
# qqplot(log2.matrix['NOTCH3',h.name],log2.matrix['NOTCH3',l.sample])
# 


######## Check DE analysis of other classical notch target genes #########

# HES1 <- 'ENSG00000114315'
# HES2 <- 'ENSG00000069812'
# HES3 <- 'ENSG00000173673'
# HES5 <- 'ENSG00000197921'
# HES6 <- 'ENSG00000144485'
# HES7 <- 'ENSG00000179111'
# 
# HEY1 <- 'ENSG00000164683'
# HEY2 <- 'ENSG00000135547'
# HEYL <- 'ENSG00000163909'
# 
# 
# load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
# de.res.metastasis.liver.vs.breast.basal[HES1,]
# de.res.metastasis.liver.vs.breast.basal[HES2,]
# 
# de.res.metastasis.liver.vs.breast.basal[HES3,]
# de.res.metastasis.liver.vs.breast.basal[HES5,]
# de.res.metastasis.liver.vs.breast.basal[HES6,]
# 
# de.res.metastasis.liver.vs.breast.basal[HES7,]
# de.res.metastasis.liver.vs.breast.basal[HEY1,]
# de.res.metastasis.liver.vs.breast.basal[HEY2,]
# de.res.metastasis.liver.vs.breast.basal[HEYL,]
# de.res.metastasis.liver.vs.breast.basal[HEY1,]
# de.res.metastasis.liver.vs.breast.basal[HEY2,]






######################### Trash code ####################################
# ############# analyze TCGA mutation data ############
# load('~/Project/InSilicoCRISPR/client-side/output/organize.TCGA.mutation.data.R.output/organize.TCGA.mutation.data.RData')
# 
# BRCA.mutation <- TCGA.mutation.data.list$BRCA
# 
# HES4.expr     <- log2.tpm.matrix[HES4,pure.TCGA.breast.cancer.polyA.Basal.sample] %>% sort
# l.sample      <- names(HES4.expr)[HES4.expr <  3]
# h.sample      <- names(HES4.expr)[HES4.expr >=  3]
# l.sample      <- intersect(BRCA.mutation %>% rownames,l.sample)
# h.sample      <- intersect(BRCA.mutation %>% rownames,h.sample)
# 
# f1            <- apply(BRCA.mutation[l.sample,],2,sum) / length(l.sample)
# f2            <- apply(BRCA.mutation[h.sample,],2,sum) / length(h.sample)
# 
# diff          <- abs(f1 - f2) %>% sort(decreasing = TRUE)
# 
# 
# 
# df.1 <- data.frame(HES4=log2.tpm.matrix[HES4,l.sample],NOTCH3=log2.tpm.matrix[NOTCH3,l.sample],m=BRCA.mutation[l.sample,'NOTCH3'])
# df.2 <- data.frame(HES4=log2.tpm.matrix[HES4,h.sample],NOTCH3=log2.tpm.matrix[NOTCH3,h.sample],m=BRCA.mutation[h.sample,'NOTCH3'])
# ggplot(rbind(df.1,df.2),aes(x=NOTCH3,y=HES4,color=factor(m))) + geom_point()


# rep.HSE4.h.up.gene    <- rownames(res)[res$log2FoldChange >  1 & res$padj < 0.05]
# rep.HSE4.h.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = rep.HSE4.h.up.gene)
# rep.HES4.h.dn.gene    <- rownames(res)[res$log2FoldChange < - 1 & res$padj < 0.05]
# rep.HES4.h.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = rep.HES4.h.dn.gene)
# 
# 
# 
# #WNT5B governs the phenotype of basal-like breast cancer by activating WNT signaling
# c.HES4.h.up.gene <- intersect(HES4.h.up.gene.df$SYMBOL,  rep.HSE4.h.up.gene.df$SYMBOL) 
# c.HES4.h.dn.gene <- intersect(HES4.h.dn.gene.df$SYMBOL,  rep.HES4.h.dn.gene.df$SYMBOL) 
# 




# basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
# basal.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = basal.up.gene)
# basal.dn.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.dn.csv", stringsAsFactors=FALSE)$x
# basal.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = basal.dn.gene)
# 
# intersect(c.HES4.h.up.gene,basal.up.gene.df$SYMBOL)
# intersect(c.HES4.h.dn.gene,basal.dn.gene.df$SYMBOL)


# ######### analyze  rppa data #########
# require(data.table)
# data <- fread('client-side/Data/brca_tcga/data_rppa_Zscores.txt',header = TRUE)
# row.names <- data$Composite.Element.REF
# data$Composite.Element.REF <- NULL
# rppa.matrix <- as.matrix(data)
# 
# rownames(rppa.matrix) <- row.names
# l.rppa.matrix <- rppa.matrix[,intersect(l.sample,colnames(rppa.matrix))]
# h.rppa.matrix <- rppa.matrix[,intersect(h.sample,colnames(rppa.matrix))]
# 
# p.value.vec <- foreach(g=rownames(l.rppa.matrix),.combine='c') %do% {
#     wilcox.test(l.rppa.matrix[g,],h.rppa.matrix[g,])$p.value  
# }
# names(p.value.vec) <- rownames(l.rppa.matrix)
# 
# 


# ######### analyze CNV data, compare CNV between HES4-high and HES4-low, TCGA #########
# HES4.expr     <- log2.tpm.matrix[HES4,pure.TCGA.breast.cancer.polyA.Basal.sample] %>% sort
# l.sample      <- names(HES4.expr)[HES4.expr <  3]
# h.sample      <- names(HES4.expr)[HES4.expr >=  3]
# 
# require(data.table)
# data                    <- fread('client-side/Data/brca_tcga/data_CNA.txt',header = TRUE)
# row.names               <- data$Hugo_Symbol
# data$Entrez_Gene_Id     <- NULL
# data$Hugo_Symbol        <- NULL
# GISTIC.matrix           <- as.matrix(data)
# rownames(GISTIC.matrix) <- row.names
# GISTIC.matrix           <- GISTIC.matrix[,colnames(GISTIC.matrix) %in% c(l.sample,h.sample)]
# 
# 
# amp.sample <- colnames(GISTIC.matrix)[GISTIC.matrix['NOTCH3',] == 2]
# l.sample   <- intersect(l.sample,colnames(GISTIC.matrix))
# h.sample   <- intersect(h.sample,colnames(GISTIC.matrix))
# 
# notch3.er.p.value <- 1- phyper(intersect(h.sample,amp.sample) %>% length(),
#        length(amp.sample),
#        length(l.sample) + length(h.sample) - length(amp.sample),
#        length(h.sample)
#        )
# 
# df.1 <- data.frame(HES4=log2.tpm.matrix[HES4,l.sample],NOTCH3=log2.tpm.matrix[NOTCH3,l.sample],m=GISTIC.matrix['NOTCH3',l.sample],n=GISTIC.matrix['NOTCH1',l.sample],o=GISTIC.matrix['NOTCH4',l.sample])
# df.2 <- data.frame(HES4=log2.tpm.matrix[HES4,h.sample],NOTCH3=log2.tpm.matrix[NOTCH3,h.sample],m=GISTIC.matrix['NOTCH3',h.sample],n=GISTIC.matrix['NOTCH1',h.sample],o=GISTIC.matrix['NOTCH4',h.sample])
# ggplot(rbind(df.1,df.2),aes(x=NOTCH3,y=HES4,color=ifelse(o == 2,'amp','no') %>% factor)) + geom_point(size=4)
# dd <- rbind(df.1,df.2)
# ggplot(dd,aes(y=HES4,x=ifelse(m == 2,'amp','no') %>% factor)) + geom_boxplot(size=2) + geom_point(size=2)
# 
# 
# data                      <- fread('client-side/Data/brca_tcga/data_linear_CNA.txt',header = TRUE)
# row.names                 <- data$Hugo_Symbol
# data$Entrez_Gene_Id       <- NULL
# data$Hugo_Symbol          <- NULL
# TCGA.cnv.matrix           <- as.matrix(data)
# rownames(TCGA.cnv.matrix) <- row.names
# TCGA.cnv.matrix           <- TCGA.cnv.matrix[,colnames(TCGA.cnv.matrix) %in% c(l.sample,h.sample)]
# 
# 
# df <- foreach(g=rownames(TCGA.cnv.matrix), .combine='rbind') %do%  {
#     p.value <- wilcox.test(TCGA.cnv.matrix[g,l.sample],TCGA.cnv.matrix[g,h.sample])$p.value  
#     ef.size <- median(TCGA.cnv.matrix[g,l.sample]) - median(TCGA.cnv.matrix[g,h.sample])
#     data.frame(ef.size=ef.size,p.value=p.value)
# }
# rownames(df) <- rownames(TCGA.cnv.matrix)
# 
# 
# ############# analyze  MET500 data #########
# load('~/Project/Cancer2CellLine/client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.cnv.RData')
# load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
# load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
# require(CNTools)
# 
# 
# 
# ######     Let us first solve the gene symbol issues. It is really nasty! ############
# ########## Finally, map all CCLE captured gene symbols to the newest HGNC version   ################
# hgnc.complete.set               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/hgnc_complete_set.txt", stringsAsFactors=FALSE)
# hgnc.complete.set               <- subset( hgnc.complete.set,select=c('hgnc_id','symbol'))
# colnames(hgnc.complete.set)     <- c('HGNC.ID','Symbol')
# df                              <- read.delim("~/Project/Cancer2CellLine/client-side/Data/somatic.mutation/CCLE.captured.gene.txt", stringsAsFactors=FALSE)
# df                              <- df[df$HGNC.Symbol != '',]
# idx                             <- match(x = df$HGNC.ID,table = hgnc.complete.set$HGNC.ID)
# mapping.df                      <- data.frame(previous.hgnc.symbol=df$HGNC.Symbol,new.hgnc.symbol=hgnc.complete.set$Symbol[idx]) %>% unique # There are duplicates in the CCLE.captured.gene.txt: two genes correspond to HGNC symbol MECOM
# rownames(mapping.df)            <- mapping.df$previous.hgnc.symbol %>% as.character
# mapping.df$previous.hgnc.symbol <- as.character(mapping.df$previous.hgnc.symbol)
# mapping.df$new.hgnc.symbol      <- as.character(mapping.df$new.hgnc.symbol) # The mapping.df is used to map all HGNC symbols to newest version
# 
# 
# hg19.RefSeq.gene.coordinates <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/hg19.RefSeq.gene.coordinates.txt", stringsAsFactors=FALSE)
# flag                         <- grepl(x=hg19.RefSeq.gene.coordinates$chrom,pattern = '_')
# hg19.RefSeq.gene.coordinates <- hg19.RefSeq.gene.coordinates[!flag,] #Let us remove the genes which are on the un-assembled contigs!
# tmp                          <- ddply(hg19.RefSeq.gene.coordinates,.(name2),function(x) data.frame(chrom=x$chrom[1], start=min(x$txStart),end=max(x$txEnd),geneid=x$name2[1],genename=x$name2[1]))
# tmp$name2                    <- NULL
# hg19.gene.info               <- tmp
# hg19.gene.info$chrom         <- gsub(x=hg19.gene.info$chrom,pattern='chr',replacement = '')
# hg19.gene.info$chrom         <- as.character(hg19.gene.info$chrom)
# hg19.gene.info$geneid        <- as.character(hg19.gene.info$geneid)
# hg19.gene.info$genename      <- as.character(hg19.gene.info$genename)
# flag                         <- which(hg19.gene.info$genename %in% mapping.df$previous.hgnc.symbol)
# for(i in flag){
#   pre.symbol <- hg19.gene.info$genename[i]  %>% as.character
#   new.symbol <- mapping.df[pre.symbol,'new.hgnc.symbol']
#   hg19.gene.info[i,'genename']        <- new.symbol
#   hg19.gene.info[i,'geneid']          <- new.symbol
# } #In hg19 RefSeq release, some CCLE captured genes are still with expired HGNC symbol, just replace them
# 
# get.id <- function(x){
#     l <- strsplit(x = x,split='\\.')  %>% unlist 
#     l[1]
# }
# flag                           <- grepl(x=MET500.sample.meta$cancer.type,pattern='Breast Invasive Ductal Carcinoma') 
# breast.cancer.sample.MET500.id <- MET500.sample.meta$MET500.id[flag] %>% as.character %>% unique
# MET500.cnv.data               <- read.csv("~/Project/Cancer2CellLine/client-side/Data/CNV/cnv_v4.csv", stringsAsFactors=FALSE)
# MET500.cnv.data               <- MET500.cnv.data[,c('Pipeline_ID','Chr','Start','End','Log2_Coverage_Ratio')]
# colnames(MET500.cnv.data)     <- c('ID','chrom','loc.start','loc.end','seg.mean')
# MET500.cnv.data$ID            <- sapply(MET500.cnv.data$ID,get.id)
# MET500.breast.cancer.cnv.data <- MET500.cnv.data[MET500.cnv.data$ID %in% breast.cancer.sample.MET500.id,]
# MET500.breast.cancer.seg      <- CNSeg(MET500.breast.cancer.cnv.data)
# MET500.breast.cancer.gene.cnv <- getRS(object = MET500.breast.cancer.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "max")@rs
# MET500.breast.cancer.gene.cnv <- MET500.breast.cancer.gene.cnv[,5:ncol(MET500.breast.cancer.gene.cnv)]
# rownames(MET500.breast.cancer.gene.cnv) <- MET500.breast.cancer.gene.cnv$genename
# MET500.breast.cancer.gene.cnv$genename  <- NULL
# MET500.breast.cancer.gene.cnv           <- as.matrix(MET500.breast.cancer.gene.cnv)
# 
# 
# TCGA.cnv.data               <- read.table("client-side/Data/gdac.broadinstitute.org_BRCA.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", stringsAsFactors=FALSE,header = TRUE)
# TCGA.cnv.data               <- TCGA.cnv.data[,c('Sample','Chromosome','Start','End','Segment_Mean')]
# colnames(TCGA.cnv.data)     <- c('ID','chrom','loc.start','loc.end','seg.mean')
# 
# get.TCGA.sample.id <- function(x) {
#     tmp <- strsplit(x=x,split='-') %>% unlist 
#     s   <-paste(tmp[1:4],collapse = '-',sep='')
#     s   <- gsub(x = s,pattern = 'A$',replacement = '')
#     s   <- gsub(x = s,pattern = 'B$',replacement = '')
#     s
# }
# TCGA.cnv.data$ID <- sapply(TCGA.cnv.data$ID,get.TCGA.sample.id)
# TCGA.seg         <- CNSeg(TCGA.cnv.data)
# TCGA.breast.cancer.gene.cnv <- getRS(object = TCGA.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "max")@rs
# TCGA.breast.cancer.gene.cnv <- TCGA.breast.cancer.gene.cnv[,5:ncol(TCGA.breast.cancer.gene.cnv)]
# rownames(TCGA.breast.cancer.gene.cnv) <- TCGA.breast.cancer.gene.cnv$genename
# TCGA.breast.cancer.gene.cnv$genename  <- NULL
# TCGA.breast.cancer.gene.cnv           <- as.matrix(TCGA.breast.cancer.gene.cnv)
# 
# 
# 
# 
# x <- MET500.sample.meta[MET500.breast.cancer.polyA.Basal.sample,]
# x <- x[x$biopsy.site == 'LIVER',]
# met500.id <- x$MET500.id %>% as.character()
# 
# MET500.breast.cancer.gene.cnv['NOTCH3',met500.id]
# 
# boxplot(MET500.breast.cancer.gene.cnv['NOTCH3',met500.id],TCGA.breast.cancer.gene.cnv['NOTCH3',h.sample],TCGA.cnv.matrix['NOTCH3',l.sample])
# 
# boxplot(MET500.breast.cancer.gene.cnv['NOTCH3',met500.id],TCGA.breast.cancer.gene.cnv['NOTCH3',c(h.sample,l.sample)])
# 
# 
# rs <- foreach(g=intersect(rownames(TCGA.cnv.matrix),rownames(MET500.breast.cancer.gene.cnv)),.combine='rbind') %do% {
#     p.value <-   wilcox.test(MET500.breast.cancer.gene.cnv[g,met500.id],TCGA.breast.cancer.gene.cnv[g,c(h.sample,l.sample)])$p.value
#     ef.size <-   median(MET500.breast.cancer.gene.cnv[g,met500.id]) - median(TCGA.breast.cancer.gene.cnv[g,c(h.sample,l.sample)])
#     data.frame(ef.size=ef.size,p.value=p.value)
# }
# rownames(rs) <- intersect(rownames(TCGA.cnv.matrix),rownames(MET500.breast.cancer.gene.cnv))


# ######### analyze methylation data #########
# 
# data                 <- fread('client-side/Data/brca_tcga/data_methylation_hm450.txt',header = TRUE)
# row.names            <- data$Hugo_Symbol
# data$Entrez_Gene_Id  <- NULL
# data$Hugo_Symbol     <- NULL
# met.matrix           <- as.matrix(data)
# rownames(met.matrix) <- row.names
# 
# l.met.matrix <- met.matrix[,intersect(l.sample,colnames(met.matrix))]
# h.met.matrix <- met.matrix[,intersect(h.sample,colnames(met.matrix))]
# 
# rs.df <- foreach(g=rownames(met.matrix),.combine='rbind') %do% {
#   p.value <- wilcox.test(l.met.matrix[g,],h.met.matrix[g,])$p.value  
#   e.size  <- median(h.met.matrix[g,]) - median(l.met.matrix[g,])
#   data.frame(effect.size=e.size,p.value=p.value)
# }
# rownames(rs.df) <- rownames(met.matrix)






########## LumB subtype ###########
# HES4 <- 'ENSG00000188290'
# 
# HES4.expr <- log2.tpm.matrix[HES4,pure.TCGA.breast.cancer.polyA.LumB.sample] %>% sort
# hist(HES4.expr,breaks=20)
# 
# l.sample <- names(HES4.expr)[HES4.expr <  2.5]
# h.sample <- names(HES4.expr)[HES4.expr >=  2.5]
# l.expr.matrix <- log2.tpm.matrix[,l.sample]
# h.expr.matrix <- log2.tpm.matrix[,h.sample]
# 
# protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
# g1                     <- get.expressed.gene(l.expr.matrix)
# g2                     <- get.expressed.gene(h.expr.matrix)
# expressed.gene         <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
# 
# expr.matrix    <- cbind(log2.read.count.matrix[expressed.gene,h.sample],log2.read.count.matrix[expressed.gene,l.sample])
# expr.matrix    <- 2^expr.matrix - 1
# 
# purity.h       <- tumor.purity.based.on.cell.line.vec[h.sample]
# purity.l       <- tumor.purity.based.on.cell.line.vec[l.sample]
# 
# df             <- data.frame(condition=c(rep(x='HES4.high',times=length(h.sample)), rep(x='HES4.low',times=length(l.sample))))
# df$purity      <- c(purity.h,purity.l)
# df$condition   <- factor(df$condition,levels = c('HES4.low','HES4.high'))
# 
# dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                          colData = df,
#                                          design = ~ condition + purity)
# dds <- DESeq(dds)
# res <- results(dds,contrast = c('condition','HES4.high','HES4.low')) %>% as.data.frame
# res <- res[order(res$pvalue),]
# res <- res[complete.cases(res),]     
# res <- as.data.frame(res)
# 
# 
# up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.05]
# up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)
# 
# dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]
# dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = dn.gene)


########## HER2 subtype ########### hmm, seems to up-regulate liver specific genes
# HES4 <- 'ENSG00000188290'
# 
# HES4.expr <- log2.tpm.matrix[HES4,pure.TCGA.breast.cancer.polyA.Her2.sample] %>% sort
# hist(HES4.expr,breaks=20)
# quantile(HES4.expr,probs = seq(0, 1, 0.05))  %>% plot
# 
# l.sample <- names(HES4.expr)[HES4.expr <  3]
# h.sample <- names(HES4.expr)[HES4.expr >=  3]
# l.expr.matrix <- log2.tpm.matrix[,l.sample]
# h.expr.matrix <- log2.tpm.matrix[,h.sample]
# 
# protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
# g1                     <- get.expressed.gene(l.expr.matrix)
# g2                     <- get.expressed.gene(h.expr.matrix)
# expressed.gene         <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
# 
# expr.matrix    <- cbind(log2.read.count.matrix[expressed.gene,h.sample],log2.read.count.matrix[expressed.gene,l.sample])
# expr.matrix    <- 2^expr.matrix - 1
# 
# purity.h       <- tumor.purity.based.on.cell.line.vec[h.sample]
# purity.l       <- tumor.purity.based.on.cell.line.vec[l.sample]
# 
# df             <- data.frame(condition=c(rep(x='HES4.high',times=length(h.sample)), rep(x='HES4.low',times=length(l.sample))))
# df$purity      <- c(purity.h,purity.l)
# df$condition   <- factor(df$condition,levels = c('HES4.low','HES4.high'))
# 
# dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                          colData = df,
#                                          design = ~ condition + purity)
# dds <- DESeq(dds)
# res <- results(dds,contrast = c('condition','HES4.high','HES4.low')) %>% as.data.frame
# res <- res[order(res$pvalue),]
# res <- res[complete.cases(res),]     
# res <- as.data.frame(res)
# 
# 
# up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.05]
# up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)
# 
# dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]
# dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = dn.gene)



# ########## LumA subtype ###########
# HES4 <- 'ENSG00000188290'
# 
# HES4.expr <- log2.tpm.matrix[HES4,pure.TCGA.breast.cancer.polyA.LumA.sample] %>% sort
# hist(HES4.expr,breaks=40)
# 
# l.sample <- names(HES4.expr)[HES4.expr <  3]
# h.sample <- names(HES4.expr)[HES4.expr >=  3]
# l.expr.matrix <- log2.tpm.matrix[,l.sample]
# h.expr.matrix <- log2.tpm.matrix[,h.sample]
# 
# protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
# g1                     <- get.expressed.gene(l.expr.matrix)
# g2                     <- get.expressed.gene(h.expr.matrix)
# expressed.gene         <- intersect(protein.coding.gene.id,c(g1,g2) %>% unique)
# 
# expr.matrix    <- cbind(log2.read.count.matrix[expressed.gene,h.sample],log2.read.count.matrix[expressed.gene,l.sample])
# expr.matrix    <- 2^expr.matrix - 1
# 
# purity.h       <- tumor.purity.based.on.cell.line.vec[h.sample]
# purity.l       <- tumor.purity.based.on.cell.line.vec[l.sample]
# 
# df             <- data.frame(condition=c(rep(x='HES4.high',times=length(h.sample)), rep(x='HES4.low',times=length(l.sample))))
# df$purity      <- c(purity.h,purity.l)
# df$condition   <- factor(df$condition,levels = c('HES4.low','HES4.high'))
# 
# dds            <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                          colData = df,
#                                          design = ~ condition + purity)
# dds <- DESeq(dds)
# res <- results(dds,contrast = c('condition','HES4.high','HES4.low')) %>% as.data.frame
# res <- res[order(res$pvalue),]
# res <- res[complete.cases(res),]     
# res <- as.data.frame(res)
# 
# 
# up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.05]
# up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.gene)
# 
# dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.05]
# dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = dn.gene)
# 


# ########## run the ssGSEA analysis ###########
# 
# require(GSVA)
# require(GSA)
# require(org.Hs.eg.db)
# 
# common.genes                     <- expressed.gene
# ensemble.to.entrez.mapping       <- revmap(org.Hs.egENSEMBL) %>% as.list
# common.genes                     <- common.genes[common.genes %in% names(ensemble.to.entrez.mapping)]
# gene.id.list                     <- ensemble.to.entrez.mapping[common.genes]
# l                                <- sapply(gene.id.list,length)
# common.genes                     <- common.genes[l == 1]
# combined.expr.matrix             <- log2.tpm.matrix[common.genes,c(l.sample,h.sample)]
# rownames(combined.expr.matrix)   <- ensemble.to.entrez.mapping[common.genes]
# 
# 
# msigdb          <-  GSA.read.gmt("~/Project/Cancer2CellLine/client-side/meta.data/c6.all.v6.1.entrez.gmt")
# genesets        <-  msigdb$genesets
# names(genesets) <-  msigdb$geneset.names 
# oncogenic.geneset.gsea.results    <-  gsva(combined.expr.matrix, genesets, method = 'ssgsea',ssgsea.norm=FALSE) #ggsea
# 
# 
# msigdb          <-  GSA.read.gmt("~/Project/Cancer2CellLine/client-side/meta.data/h.all.v6.1.entrez.gmt")
# genesets        <-  msigdb$genesets
# names(genesets) <-  msigdb$geneset.names 
# hallmark.geneset.gsea.results    <-  gsva(combined.expr.matrix , genesets, method = 'ssgsea',ssgsea.norm=FALSE) #ggsea
# 
# rs.df <- foreach(gs = rownames(hallmark.geneset.gsea.results),.combine = 'rbind') %do% {
#     p.value <- wilcox.test(hallmark.geneset.gsea.results[gs,h.sample],hallmark.geneset.gsea.results[gs,l.sample])$p.value
#     e.size <- median(hallmark.geneset.gsea.results[gs,h.sample]) - median(hallmark.geneset.gsea.results[gs,l.sample])
#     data.frame(effect.size=e.size,p.value=p.value)
# }
# rownames(rs.df) <- rownames(hallmark.geneset.gsea.results)
# rs.df$fdr <- p.adjust(rs.df$p.value,method='fdr')
# 
# 
# rs.df <- foreach(gs = rownames(oncogenic.geneset.gsea.results),.combine = 'rbind') %do% {
#   p.value <- wilcox.test(oncogenic.geneset.gsea.results[gs,h.sample],oncogenic.geneset.gsea.results[gs,l.sample])$p.value
#   e.size <- median(oncogenic.geneset.gsea.results[gs,h.sample]) - median(oncogenic.geneset.gsea.results[gs,l.sample])
#   data.frame(effect.size=e.size,p.value=p.value)
# }
# rownames(rs.df) <- rownames(oncogenic.geneset.gsea.results)
# rs.df$fdr <- p.adjust(rs.df$p.value,method='fdr')
####### compare the TME ##########



####### compare the mutation ########



####### compare the CNV ########



######### compare methylation ##########



# rs <- foreach(g= g.vec,.combine='rbind') %do% {
#     p.value <- wilcox.test(l.expr.matrix[g,],h.expr.matrix[g,])$p.value
#     effect.size <- median(h.expr.matrix[g,]) - median(l.expr.matrix[g,])
#     data.frame(effect.size=effect.size,p.value=p.value)
# }
# rownames(rs) <- g.vec
# rs <- rs[order(rs$p.value),]
# rs$fdr <- p.adjust(rs$p.value,method='fdr')
# rs <- rs[complete.cases(rs),]
# 
# 
# 
# 
# 
# 
# 
# load('~/Project/GeneFishing/RData/normalized.GTex.with.ncRNA.median.larger.than.0.1.RData')
# HES4 <- 'ENSG00000188290'

#CEBPD <- 'ENSG00000221869'
#HES4 <- 'ENSG00000221869' # CEBPD
#HES4 <- 'ENSG00000164916'# FOXK1
#HES4 <- 'ENSG00000008128' #CDK11A
#HES4 <- 'ENSG00000164916' # FOXK1
#HES4 <- 'ENSG00000104885' #DOT1L
#HES4 <- 'ENSG00000167491' #GATAD2A

############# breast cancer ############# 
# data           <- normalized.GTex.data[['Breast - Mammary Tissue']]
# cor.matrix     <- cor(x=data,y=data[,HES4],method='spearman')
# cor.vec        <- c(cor.matrix)
# names(cor.vec) <- rownames(cor.matrix)
# cor.vec        <- sort(cor.vec,decreasing = TRUE)
# up.outlier     <- quantile(cor.vec)[4] + 1.5 * IQR(cor.vec)
# dn.outlier     <- quantile(cor.vec)[1] - 1.5 * IQR(cor.vec)
# 
# target.gene       <- names(cor.vec)[cor.vec > up.outlier]
# target.gene.df    <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = target.gene)
# br.target.gene.df <- target.gene.df[complete.cases(target.gene.df),]
# 
# luma.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/luma.up.csv",  stringsAsFactors=FALSE)$x
# lumb.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv",  stringsAsFactors=FALSE)$x
# basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
# her2.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv",  stringsAsFactors=FALSE)$x
# 
# luma.up.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = luma.up.gene)
# lumb.up.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = lumb.up.gene)
# basal.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = basal.up.gene)
# her2.up.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = her2.up.gene)
# 
# 
# intersect(br.target.gene.df$SYMBOL,luma.up.gene.df$SYMBOL)
# intersect(br.target.gene.df$SYMBOL,lumb.up.gene.df$SYMBOL)
# intersect(br.target.gene.df$SYMBOL,her2.up.gene.df$SYMBOL)
# intersect(br.target.gene.df$SYMBOL,basal.up.gene.df$SYMBOL)
# 
# 
# ### genes up-regulated in at least two DE comparisions
# #tmp         <- c(luma.up.gene,lumb.up.gene) %>% unique
# tmp         <- c(luma.up.gene,lumb.up.gene,basal.up.gene,her2.up.gene)
# freq.table  <- table(tmp) %>% as.data.frame
# ov.gene     <- freq.table$tmp[freq.table$Freq >= 2] %>% as.character()
# ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = ov.gene)
# ov.gene.df  <- ov.gene.df[complete.cases(ov.gene.df),]
# intersect(ov.gene.df$SYMBOL,target.gene.df$SYMBOL)
# 
# ### genes up-regulated in all four DE comparisions
# tmp              <- c(luma.up.gene,lumb.up.gene,basal.up.gene,her2.up.gene)
# freq.table       <- table(tmp) %>% as.data.frame
# high.ov.gene     <- freq.table$tmp[freq.table$Freq >= 4] %>% as.character()
# high.ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = high.ov.gene)
# high.ov.gene.df  <- high.ov.gene.df[complete.cases(high.ov.gene.df),]
# intersect(high.ov.gene.df$SYMBOL,target.gene.df$SYMBOL)
# 
# 
# ################## prostate cancer ##############
# pr.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.prostate.cancer.R.output/up.gene.csv", stringsAsFactors=FALSE)$x
# pr.up.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = pr.up.gene)
# data           <- normalized.GTex.data[['Prostate']]
# cor.matrix     <- cor(x=data,y=data[,HES4],method='spearman')
# cor.vec        <- c(cor.matrix)
# names(cor.vec) <- rownames(cor.matrix)
# cor.vec        <- sort(cor.vec,decreasing = TRUE)
# up.outlier     <- quantile(cor.vec)[4] + 1.5 * IQR(cor.vec)
# 
# target.gene       <- names(cor.vec)[cor.vec > up.outlier]
# target.gene.df    <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = target.gene)
# pr.target.gene.df <- target.gene.df[complete.cases(target.gene.df),]
# 
# 
# 
# intersect(pr.up.gene.df$SYMBOL,pr.target.gene.df$SYMBOL)



# #################################
# load('server-side/RData/Breast Invasive Carcinoma.RData')
# load('client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')
# 
# s <- intersect(names(tumor.purity.based.on.cell.line.vec),colnames(log2.tpm.matrix))
# plot(x=tumor.purity.based.on.cell.line.vec[s],y=log2.tpm.matrix[HES4,s])
# 
# 
# cor.matrix     <- cor(x=log2.tpm.matrix[HES4,],y=log2.tpm.matrix %>% t,method='spearman')
# cor.vec        <- c(cor.matrix)
# names(cor.vec) <- colnames(cor.matrix)
# cor.vec        <- sort(cor.vec,decreasing = TRUE)
# 
# 
# HES4.expr <- log2.tpm.matrix[HES4,] %>% sort
# 
# l.sample <- names(HES4.expr)[1:100]
# h.sample <- names(HES4.expr)[992:1092]
# l.expr.matrix <- log2.tpm.matrix[,l.sample]
# h.expr.matrix <- log2.tpm.matrix[,h.sample]
# 
# protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
# g.vec                  <- intersect(protein.coding.gene.id,rownames(l.expr.matrix))
# rs <- foreach(g= g.vec,.combine='rbind') %do% {
#     p.value <- wilcox.test(l.expr.matrix[g,],h.expr.matrix[g,])$p.value  
#     effect.size <- median(h.expr.matrix[g,]) - median(l.expr.matrix[g,])
#     data.frame(effect.size=effect.size,p.value=p.value)
# }
# rownames(rs) <- g.vec
# rs <- rs[order(rs$p.value),]
# rs$fdr <- p.adjust(rs$p.value,method='fdr')
# rs <- rs[complete.cases(rs),]



# ###################3
# require(affy) # GSE5460 GSE3744
# setwd('~/Project/BreastCancerMetaPotenial/client-side/Data/GSE5460_RAW/')
# cel.file   <- c(system('ls | grep cel.gz',intern = TRUE,wait = TRUE),system('ls | grep CEL.gz',intern = TRUE,wait = TRUE)) %>% unique
# cel.file.1 <- sapply(cel.file,function(x) paste('/Users/liuke/Project/BreastCancerMetaPotenial/client-side/Data/GSE5460_RAW/',x,sep=''))
# 
# setwd('~/Project/BreastCancerMetaPotenial/client-side/Data/GSE3744_RAW/')
# cel.file   <- c(system('ls | grep cel.gz',intern = TRUE,wait = TRUE),system('ls | grep CEL.gz',intern = TRUE,wait = TRUE)) %>% unique
# cel.file.2 <- sapply(cel.file,function(x) paste('/Users/liuke/Project/BreastCancerMetaPotenial/client-side/Data/GSE3744_RAW/',x,sep=''))
# 
# cel.file   <- c(cel.file.1)
# 
# raw.data                      <- ReadAffy(filenames=cel.file)
# data.rma.norm                 <- rma(raw.data)
# y                             <- exprs(data.rma.norm)
# x                             <- colnames(y)
# pattern                       <- 'GSM\\d+'
# m                             <- regexpr(pattern, x)
# colnames(y)                   <- regmatches(x, m)
# microarray.breast.cancer.data <- y
# 
# HES4.probe.set     <- '227347_x_at'
# tmp                <- cor(microarray.breast.cancer.data[HES4.probe.set,], microarray.breast.cancer.data %>% t, method='spearman')
# microarray.cor.vec <- c(tmp)
# names(microarray.cor.vec) <- colnames(tmp)
# microarray.cor.vec        <- sort(microarray.cor.vec,decreasing = TRUE)
# 
# require(hgu133plus2.db)
# require(AnnotationDbi)
# 
# keytypes(hgu133plus2.db)
# 
# probe.to.ensemble.mapping.df <- AnnotationDbi::select(x=hgu133plus2.db,keys = names(microarray.cor.vec),keytype = 'PROBEID',columns = 'ENSEMBL')
# probe.to.ensemble.mapping.df <- probe.to.ensemble.mapping.df[complete.cases(probe.to.ensemble.mapping.df),]
# 
# mapping.cnt                  <- ddply(probe.to.ensemble.mapping.df,.(PROBEID),nrow)
# unique.probe                 <- mapping.cnt$PROBEID[mapping.cnt$V1 == 1]
# probe.to.ensemble.mapping.df <- probe.to.ensemble.mapping.df[probe.to.ensemble.mapping.df$PROBEID %in% unique.probe ,]
# 
# microarray.df                 <- cbind(probe.to.ensemble.mapping.df,array.cor.value=microarray.cor.vec[probe.to.ensemble.mapping.df$PROBEID %>% as.character()])
# tmp                           <- ddply(microarray.df,.(ENSEMBL),function(x) mean(x$array.cor.value))
# aggregated.microarray.cor.vec <- tmp$V1
# names(aggregated.microarray.cor.vec) <- tmp$ENSEMBL
# 
# setwd('/Users/liuke/Project/BreastCancerMetaPotenial')
# # # require(sva)
# # # 
# # # batch.vec <- ifelse(grepl(x=colnames(y),pattern='GSM85') ,1,2) 
# # # combat_edata1 = ComBat(dat=y, batch=batch.vec, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
# # HES4.probe.set <- '227347_x_at'
# # hist(y[HES4.probe.set,],breaks=20)
# # 
# # combat_edata1[HES4.probe.set,] %>% sort %>% plot
# # 
# # combat_edata1 <- combat_edata1[,order(combat_edata1[HES4.probe.set,])]
# # combat_edata1 <- combat_edata1[,c(1:37,137:174)]
# # 
# # col.data           <- data.frame(HES4.status = ifelse(combat_edata1[HES4.probe.set,] > 8.25,'HES4.high','HES4.low'))
# # rownames(col.data) <- colnames(combat_edata1)
# # col.data$HES4.status    <- factor(col.data$HES4.status,levels = c('HES4.low','HES4.high'))
# # 
# # ds  <- model.matrix(~HES4.status,col.data)
# # fit <- lmFit(combat_edata1, ds)
# # fit <- eBayes(fit)
# # array.rs  <- topTable(fit, coef="HES4.statusHES4.high",number=nrow(y))
# # 
# # array.rs  <- topTable(fit, coef="HES4.statusHES4.high")
# # array.rs[array.rs$adj.P.Val < 0.05,] %>% View
# # 
# # 
# # array.rs.p.value <- foreach(g=rownames(combat_edata1),.combine='c') %do% {
# #     wilcox.test(combat_edata1[g,1:37],combat_edata1[g,38:75])$p.value
# # }
# # names(array.rs.p.value) <- rownames(combat_edata1)
# # 
# # 
# # plot(x=array.rs$logFC,y=-1 * (array.rs$adj.P.Val %>% log10))
# # 
# # 
# # 
# # cor.matrix <- cor(combat_edata1[HES4.probe.set,] %>% c,combat_edata1 %>% t,method='spearman')
# # cor.vec <- c(cor.matrix)
# # names(cor.vec) <- colnames(cor.matrix)
# # cor.vec <- sort(cor.vec,decreasing = TRUE)
# 
# 
# 
# #################################
# require(DESeq2)
# require(ggplot2)
# require(dplyr)
# load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
# load('server-side/RData//Breast Invasive Carcinoma.RData')
# load('~/Project/InSilicoCRISPR/client-side/output/organize.TCGA.clinical.data.R.output/organize.TCGA.clinical.data.RData')
# 
# 
# TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix
# TCGA.breast.cancer.log2.tpm.matrix       <- log2.tpm.matrix
# TCGA.breast.cancer.log2.tpm.matrix       <- TCGA.breast.cancer.log2.tpm.matrix[,setdiff(colnames(TCGA.breast.cancer.log2.tpm.matrix), TCGA.breast.cancer.polyA.Basal.sample)]
# 
# HES4            <- 'ENSG00000188290' # wow, bimodal distribution in Basal
# tmp             <- cor(TCGA.breast.cancer.log2.tpm.matrix[HES4,], TCGA.breast.cancer.log2.tpm.matrix %>% t, method='spearman')
# rna.seq.cor.vec <- c(tmp)
# names(rna.seq.cor.vec)    <- colnames(tmp)
# rna.seq.cor.vec <- rna.seq.cor.vec[is.na(rna.seq.cor.vec) == FALSE]
# rna.seq.cor.vec           <- sort(rna.seq.cor.vec,decreasing = TRUE)
# ################
# 
# g <- intersect(names(rna.seq.cor.vec),names(aggregated.microarray.cor.vec))
# g <- setdiff(g,HES4)
# plot(x=aggregated.microarray.cor.vec[g] ,y=rna.seq.cor.vec[g] )
# 
# x <- aggregated.microarray.cor.vec[g]  %>% scale(center = TRUE,scale = TRUE) 
# y <- rna.seq.cor.vec[g]                %>% scale(center = TRUE,scale = TRUE) 
# #names(x) <- g
# #names(y) <- g
# 
# # g[x>=1.5 & y>=1.5] %>% View
# # 
# # z <- x + y
# # z <- sort(z,decreasing = TRUE)
# # View(z)
# 
# M <- y - x
# A <- y + x
# plot(x=A,y=M)
# g[A >= 4 & abs(M) <= 1] %>% View
# g[A <= -4  & abs(M) <= 1] %>% View
# plot(x=A[9000:15000,1] %>% c,y=M[9000:15000,1] %>% c)
# 
# 
# 
# # TCGA.breast.cancer.log2.tpm.matrix <- TCGA.breast.cancer.log2.tpm.matrix[,order(TCGA.breast.cancer.log2.tpm.matrix[HES4,])]
# # 
# # TCGA.breast.cancer.log2.tpm.matrix <- TCGA.breast.cancer.log2.tpm.matrix[,c(1:90,806:896)]
# # 
# # pca.rs <- prcomp(TCGA.breast.cancer.log2.tpm.matrix %>% t)
# # plot(pca.rs$x[,1:2])
# 
# 
# #Basal.log2.tpm.matrix <- TCGA.breast.cancer.log2.tpm.matrix
# 
# #Basal.log2.tpm.matrix  <-  TCGA.breast.cancer.log2.tpm.matrix[,TCGA.breast.cancer.polyA.Her2.sample]
# #Basal.read.count.matrix <-  2^TCGA.breast.cancer.log2.read.count.matrix[,TCGA.breast.cancer.polyA.Her2.sample] - 1
# 
# 
# #Basal.log2.tpm.matrix  <-  TCGA.breast.cancer.log2.tpm.matrix[,TCGA.breast.cancer.polyA.LumB.sample]
# #Basal.read.count.matrix <-  2^TCGA.breast.cancer.log2.read.count.matrix[,TCGA.breast.cancer.polyA.LumB.sample] - 1
# 
# 
# HES4           <- 'ENSG00000188290' # wow, bimodal distribution in Basal
# df             <- data.frame(HES4.status = ifelse(TCGA.breast.cancer.log2.tpm.matrix[HES4,] >=3, 'HES4.high','HES4.low'))
# df$HES4.status <- factor(df$HES4.status,levels = c('HES4.low','HES4.high'))
# flag           <- apply(TCGA.breast.cancer.log2.tpm.matrix,1,function(x) median(x) > log2(1+0.1))
# TCGA.breast.cancer.log2.tpm.matrix <- TCGA.breast.cancer.log2.tpm.matrix[flag,]
# 
# 
# res <- foreach(g= rownames(TCGA.breast.cancer.log2.tpm.matrix),.combine='rbind') %do% {
#     p.value     <- wilcox.test(TCGA.breast.cancer.log2.tpm.matrix[g,df$HES4.status == 'HES4.high'],TCGA.breast.cancer.log2.tpm.matrix[g,df$HES4.status != 'HES4.high'])$p.value 
#     effect.size <- median(TCGA.breast.cancer.log2.tpm.matrix[g,df$HES4.status == 'HES4.high']) - median(TCGA.breast.cancer.log2.tpm.matrix[g,df$HES4.status != 'HES4.high'])
#     data.frame(effect.size=effect.size,p.value=p.value)
# }
# rownames(res) <- rownames(TCGA.breast.cancer.log2.tpm.matrix)
# res$padj      <- p.adjust(res$p.value,method='fdr')
# 
# plot(x=res$effect.size,y= -1 * log10(res$padj))
# 
# up.gene <- rownames(res)[res$padj < 0.01 & res$effect.size >  1] 
# dn.gene <- rownames(res)[res$padj < 0.01 & res$effect.size < -1] 
# 
# write.csv(x=up.gene,file='~/Desktop/up.gene.csv',quote=FALSE)
# write.csv(x=dn.gene,file='~/Desktop/dn.gene.csv',quote=FALSE)
# 
# 
# 
# 
# HES4           <- 'ENSG00000188290' # wow, bimodal distribution in Basal
# df             <- data.frame(HES4.status = ifelse(Basal.log2.tpm.matrix[HES4,] >=3, 'HES4.high','HES4.low'))
# df$HES4.status <- factor(df$HES4.status,levels = c('HES4.low','HES4.high'))
# flag           <- apply(Basal.log2.tpm.matrix,1,function(x) median(x) > log2(1+0.1))
# Basal.read.count.matrix <- Basal.read.count.matrix[flag,]
# 
# 
# dds          <- DESeqDataSetFromMatrix(countData = round(Basal.read.count.matrix),
#                                        colData = df,
#                                        design = ~ HES4.status )
# 
# dds <- DESeq(dds)
# res <- results(dds,contrast = c('HES4.status','HES4.high','HES4.low')) %>% as.data.frame
# res <- res[order(res$pvalue),]
# res <- res[complete.cases(res),]
# protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
# res <- res[rownames(res) %in% protein.coding.gene.id,]
# 
# up.gene <- rownames(res)[res$padj < 0.01 & res$log2FoldChange >  1] 
# dn.gene <- rownames(res)[res$padj < 0.01 & res$log2FoldChange < -1] 
# 
# write.csv(x=up.gene,file='~/Desktop/up.gene.csv',quote=FALSE)
# write.csv(x=dn.gene,file='~/Desktop/dn.gene.csv',quote=FALSE)
# 
# # Abstract 5002: KRT13 promotes stemness and drives metastasis in breast cancer through direct interaction with plakoglobin-desmoplakin complexes regulating c-Myc signaling pathway
# df$y <- Basal.log2.tpm.matrix['ENSG00000168298',rownames(df)]
# df <- df[df$y < 2.5,]
# ggplot(df,aes(x=HES4.status,y=y)) + geom_boxplot()

#ref http://www.discoverymedicine.com/Jing-Liu/2018/05/collagen-col1a1-promotes-metastasis-of-breast-cancer-potential-therapeutic-target/
# gene SCX (LumB), required for extra celluar matrix related gene expression https://www.ncbi.nlm.nih.gov/pubmed/18802027
#gene KRT13(Basal) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2654210/
# collegan promotes metastasis. Ref: Membrane associated collagen XIII promotes cancer metastasis and enhances anoikis resistance
#ref: CCDC85B promotes non-small cell lung cancer cell proliferation and invasion.

#Rational design of anti-GITR-based combination immunotherapy (TNFRSF18 aka known as GITR)
