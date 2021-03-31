require(CNTools)



load('client-side/output/hg19.gene.info.R.output/hg19.gene.info.RData')
################# organize TCGA PRAD cnv data #####################
data                                    <- read.table("client-side/Data/cBioPortal/prad_tcga_pan_can_atlas_2018//data_cna_hg19.seg", stringsAsFactors=FALSE,header = TRUE)
TCGA.seg                                <- CNSeg(data)
TCGA.prostate.cancer.gene.cnv           <- getRS(object = TCGA.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "max")@rs
TCGA.prostate.cancer.gene.cnv           <- TCGA.prostate.cancer.gene.cnv[,5:ncol(TCGA.prostate.cancer.gene.cnv)]
rownames(TCGA.prostate.cancer.gene.cnv) <- TCGA.prostate.cancer.gene.cnv$genename
TCGA.prostate.cancer.gene.cnv$genename  <- NULL
TCGA.prostate.cancer.cnv.matrix         <- as.matrix(TCGA.prostate.cancer.gene.cnv)



load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
################# organize MET500 cnv data #####################
get.id <- function(x){
  l <- strsplit(x = x,split='\\.')  %>% unlist 
  l[1]
}
flag                                    <- grepl(x=MET500.sample.meta$cancer.type,pattern='Prostate Adenocarcinoma') 
prostate.cancer.sample.MET500.id          <- MET500.sample.meta$MET500.id[flag] %>% as.character %>% unique
MET500.cnv.data                         <- read.csv("~/Project/Cancer2CellLine/client-side/Data/CNV/cnv_v4.csv", stringsAsFactors=FALSE)
MET500.cnv.data                         <- MET500.cnv.data[,c('Pipeline_ID','Chr','Start','End','Log2_Coverage_Ratio')]
colnames(MET500.cnv.data)               <- c('ID','chrom','loc.start','loc.end','seg.mean')
MET500.cnv.data$ID                      <- sapply(MET500.cnv.data$ID,get.id)
MET500.prostate.cancer.cnv.data           <- MET500.cnv.data[MET500.cnv.data$ID %in% prostate.cancer.sample.MET500.id,]
MET500.prostate.cancer.seg                <- CNSeg(MET500.prostate.cancer.cnv.data)
MET500.prostate.cancer.gene.cnv           <- getRS(object = MET500.prostate.cancer.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "max")@rs
MET500.prostate.cancer.gene.cnv           <- MET500.prostate.cancer.gene.cnv[,5:ncol(MET500.prostate.cancer.gene.cnv)]
rownames(MET500.prostate.cancer.gene.cnv) <- MET500.prostate.cancer.gene.cnv$genename
MET500.prostate.cancer.gene.cnv$genename  <- NULL
MET500.prostate.cancer.cnv.matrix         <- as.matrix(MET500.prostate.cancer.gene.cnv)

########################################################
save(file='client-side/output/organize.TCGA.and.MET500.prostate.cancer.cnv.data.R.output/organize.TCGA.and.MET500.prostate.cancer.cnv.data.RData',list=c('TCGA.prostate.cancer.cnv.matrix','MET500.prostate.cancer.cnv.matrix'))
