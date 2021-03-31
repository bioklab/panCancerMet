require(CNTools)
require(data.table)


load('client-side/output/hg19.gene.info.R.output/hg19.gene.info.RData')
################# organize TCGA BRCA cnv data #####################
data                                  <- read.table("client-side/Data/cBioPortal/brca_tcga/data_cna_hg19.seg", stringsAsFactors=FALSE,header = TRUE)
TCGA.seg                              <- CNSeg(data)
TCGA.breast.cancer.gene.cnv           <- getRS(object = TCGA.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "max")@rs
TCGA.breast.cancer.gene.cnv           <- TCGA.breast.cancer.gene.cnv[,5:ncol(TCGA.breast.cancer.gene.cnv)]
rownames(TCGA.breast.cancer.gene.cnv) <- TCGA.breast.cancer.gene.cnv$genename
TCGA.breast.cancer.gene.cnv$genename  <- NULL
TCGA.breast.cancer.cnv.matrix           <- as.matrix(TCGA.breast.cancer.gene.cnv)



load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
################# organize MET500 cnv data #####################
get.id <- function(x){
  l <- strsplit(x = x,split='\\.')  %>% unlist 
  l[1]
}
flag                                    <- grepl(x=MET500.sample.meta$cancer.type,pattern='Breast Invasive Ductal Carcinoma') 
breast.cancer.sample.MET500.id          <- MET500.sample.meta$MET500.id[flag] %>% as.character %>% unique
MET500.cnv.data                         <- read.csv("~/Project/Cancer2CellLine/client-side/Data/CNV/cnv_v4.csv", stringsAsFactors=FALSE)
MET500.cnv.data                         <- MET500.cnv.data[,c('Pipeline_ID','Chr','Start','End','Log2_Coverage_Ratio')]
colnames(MET500.cnv.data)               <- c('ID','chrom','loc.start','loc.end','seg.mean')
MET500.cnv.data$ID                      <- sapply(MET500.cnv.data$ID,get.id)
MET500.breast.cancer.cnv.data           <- MET500.cnv.data[MET500.cnv.data$ID %in% breast.cancer.sample.MET500.id,]
MET500.breast.cancer.seg                <- CNSeg(MET500.breast.cancer.cnv.data)
MET500.breast.cancer.gene.cnv           <- getRS(object = MET500.breast.cancer.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "max")@rs
MET500.breast.cancer.gene.cnv           <- MET500.breast.cancer.gene.cnv[,5:ncol(MET500.breast.cancer.gene.cnv)]
rownames(MET500.breast.cancer.gene.cnv) <- MET500.breast.cancer.gene.cnv$genename
MET500.breast.cancer.gene.cnv$genename  <- NULL
MET500.breast.cancer.cnv.matrix         <- as.matrix(MET500.breast.cancer.gene.cnv)

########################################################
save(file='client-side/output/organize.TCGA.and.MET500.breast.cancer.cnv.data.R.output/organize.TCGA.and.MET500.breast.cancer.cnv.data.RData',list=c('TCGA.breast.cancer.cnv.matrix','MET500.breast.cancer.cnv.matrix'))
