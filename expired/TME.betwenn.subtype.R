################## Here let us compare TME difference between LumB and Basal, liver metastasis ####################

########### First, remove gene DE between LumB and Basal type ######################
####### Here we could use single-cell RNASeq data, or CCLE cell line data

library(xCell)
library(plyr)
library(dplyr)
library (VennDiagram)
library(gplots)
library(foreach)
library(org.Hs.eg.db)
library(pheatmap)
library(foreach)
library(AnnotationDbi)
source('client-side/code/util.R')
load('~/Project/Cancer2CellLine/server-side/RData/CCLE_BRACA.RData')

######## map ENSEMBLE id to gene symbol ##############
protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
ensemble.gene.id       <- intersect(rownames(GTex.liver.log2.fpkm.matrix),protein.coding.gene.id)

df    <- AnnotationDbi::select(org.Hs.eg.db,keys=ensemble.gene.id, columns=c('SYMBOL'), keytype="ENSEMBL")
df    <- df[complete.cases(df),]
tmp   <- ddply(df,.(ENSEMBL),nrow)
tmp   <- tmp[order(tmp$V1,decreasing = TRUE),]
g1    <- tmp$ENSEMBL[tmp$V1 == 1] %>% as.character()
tmp   <- ddply(df,.(SYMBOL),nrow)
tmp   <- tmp[order(tmp$V1,decreasing = TRUE),]
g2    <- tmp$SYMBOL[tmp$V1 == 1] %>% as.character()
mapping.df           <- df[(df$ENSEMBL %in% g1)  & (df$SYMBOL %in% g2),]
rownames(mapping.df) <- mapping.df$ENSEMBL

data                  <- BRACA.log2.fpkm.matrix
data                  <- data[rownames(data) %in% mapping.df$ENSEMBL,]
rownames(data)        <- mapping.df[rownames(data),'SYMBOL']
rank.data             <- apply(data,2,rank)
CCLE.data.for.xCell   <- data

krt14 <- 'ENSG00000186847'
hist(data['KRT14',])
plot(data['KRT14',] %>% sort)



df.LumB         <-  read.xls("client-side/Data/TableS3.xls", sheet = 2, header = TRUE)
df.LumB         <- df.LumB[order(df.LumB$p.value),]
df.Basal        <-  read.xls("client-side/Data/TableS3.xls", sheet = 4, header = TRUE)
df.Basal        <- df.Basal[order(df.Basal$p.value),]

Basal.cell.line <- df.Basal$cell.line[df.Basal$fdr < 0.05] %>% as.character
LumB.cell.line  <- df.LumB$cell.line[df.LumB$fdr < 0.05] %>% as.character
get.name <- function(x)
{
    tmp <- strsplit(x=x,split='_') %>% unlist
    tmp[1]
}
Basal.cell.line       <- sapply(Basal.cell.line,get.name)
LumB.cell.line        <- sapply(LumB.cell.line,get.name)
krt14.Basal.cell.line <- data['KRT14',Basal.cell.line]
krt14.LumB.cell.line  <- data['KRT14',LumB.cell.line]

Basal.cell.line <- names(krt14.Basal.cell.line)[krt14.Basal.cell.line > 4]
LumB.cell.line  <- names(krt14.LumB.cell.line)[krt14.LumB.cell.line < 4]

Basal.rank.data <- rank.data[,Basal.cell.line]
LumB.rank.data <- rank.data[,LumB.cell.line]

rank.de.df <- foreach(i = 1:nrow(Basal.rank.data),.combine = 'rbind') %do% {
  data.frame(delta = median(Basal.rank.data[i,]) - median(LumB.rank.data[i,]),p.value=wilcox.test(Basal.rank.data[i,],LumB.rank.data[i,])$p.value) 
}
rownames(rank.de.df) <- rownames(Basal.rank.data)
rank.de.df$fdr       <- p.adjust(rank.de.df$p.value,method='fdr')
flag                 <- (rank.de.df$fdr < 0.1) 
de.gene              <- rownames(rank.de.df)[flag]

xCell.corrected.gene.set <- foreach(x=xCell.data$signatures) %do% {
  setdiff(x@geneIds,de.gene) 
}
names(xCell.corrected.gene.set) <- names(xCell.data$signatures) ############# remove DE genes from xCell signatures 

L <- list()
foreach(i = names(xCell.corrected.gene.set)) %do% {
  tmp <- strsplit(x = i,split='%')[[1]] %>% unlist      
  cell.type <- tmp[1]
  L[[cell.type]] <- c(L[[cell.type]],xCell.corrected.gene.set[[i]])
  L[[cell.type]] <- unique(L[[cell.type]])
}
xCell.corrected.gene.set <- L



CCLE.cancer.ssgsea.scores <- GSVA::gsva(expr=CCLE.data.for.xCell[,c(Basal.cell.line,LumB.cell.line)],
                                        xCell.corrected.gene.set, method = "ssgsea",
                                        ssgsea.norm = FALSE)


CCLE.p.value.vec <- foreach(i= 1:nrow(CCLE.cancer.ssgsea.scores),.combine='c') %do% {
  wilcox.test(CCLE.cancer.ssgsea.scores[i,Basal.cell.line],CCLE.cancer.ssgsea.scores[i,LumB.cell.line])$p.value  
}
names(CCLE.p.value.vec) <- rownames(CCLE.cancer.ssgsea.scores)



####################################################################################################################################
############# compute correlation values with cell line, as estimate of tumor purity
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





##########################

load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')

TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix

data                  <- cbind(TCGA.breast.cancer.log2.fpkm.matrix,MET500.log2.fpkm.matrix[rownames(TCGA.breast.cancer.log2.fpkm.matrix),])
data                  <- data[rownames(data) %in% mapping.df$ENSEMBL,]
rownames(data)        <- mapping.df[rownames(data),'SYMBOL']
cancer.data.for.xCell <- data

MET500.liver.sample   <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 


MET500.sample                              <- MET500.breast.cancer.polyA.Basal.sample
MET500.polyA.Basal.liver.metastasis.sample <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample                              <- MET500.breast.cancer.polyA.LumB.sample
MET500.polyA.LumB.liver.metastasis.sample  <- intersect(MET500.sample,MET500.liver.sample)

MET500.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(MET500.polyA.Basal.liver.metastasis.sample,MET500.polyA.LumB.liver.metastasis.sample)],
                                   xCell.corrected.gene.set, method = "ssgsea",
                                   ssgsea.norm = FALSE)


MET500.p.value.vec <- foreach(i=1:nrow(MET500.ssgsea.scores),.combine='c') %do% {
    wilcox.test(MET500.ssgsea.scores[i,MET500.polyA.Basal.liver.metastasis.sample],MET500.ssgsea.scores[i,MET500.polyA.LumB.liver.metastasis.sample])$p.value  
}
names(MET500.p.value.vec) <- rownames(MET500.ssgsea.scores)




TCGA.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.breast.cancer.polyA.Basal.sample,TCGA.breast.cancer.polyA.LumB.sample)],
                                   xCell.corrected.gene.set, method = "ssgsea",
                                   ssgsea.norm = FALSE)


TCGA.df <- foreach(i=1:nrow(TCGA.ssgsea.scores),.combine='rbind') %do% {
    p.value <- wilcox.test(TCGA.ssgsea.scores[i,TCGA.breast.cancer.polyA.Basal.sample],TCGA.ssgsea.scores[i,TCGA.breast.cancer.polyA.LumB.sample])$p.value  
    delta   <- median(TCGA.ssgsea.scores[i,TCGA.breast.cancer.polyA.Basal.sample]) - median(TCGA.ssgsea.scores[i,TCGA.breast.cancer.polyA.LumB.sample])
    data.frame(delta=delta,p.value=p.value)
}
rownames(TCGA.df) <- rownames(TCGA.ssgsea.scores)






# TCGA.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.breast.cancer.polyA.Basal.sample,TCGA.breast.cancer.polyA.LumB.sample)],
#                                  xCell.corrected.gene.set, method = "ssgsea",
#                                  ssgsea.norm = FALSE)
# 


# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500 <- rs.MET500$correlation.matrix[,'EFM192A_BREAST'] # use EFM192A cell line
# purity.TCGA   <- rs.TCGA$correlation.matrix[,'EFM192A_BREAST']




