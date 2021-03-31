require(genefu)
require(data.table)
require(dplyr)

###### organize MBC meta data ######################
tmp                       <- fread(input='client-side/Data/brca_mbcproject_wagle_2017/data_clinical_sample.txt',header = TRUE,skip=4) %>% as.data.frame
tmp                       <- tmp[,c('SAMPLE_ID','PATIENT_ID','BX_LOCATION','CANCER_TYPE_DETAILED')]
#tmp                       <- tmp[tmp$CANCER_TYPE_DETAILED =='Breast Invasive Ductal Carcinoma',]
MBC.sample.meta           <- tmp[,1:3]
colnames(MBC.sample.meta) <- c('sample.id','subject.id','biopsy.site')
rownames(MBC.sample.meta) <- MBC.sample.meta$sample.id

###### PAM50 subtypeing ###########

tmp <- fread(input='client-side/Data/brca_mbcproject_wagle_2017/data_RNA_Seq_v2_expression_median.txt') %>% as.data.frame

pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') # well, I borrow some information from the metastatic breast cancer evaluation project 
pam50.gene                  <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)


pam50.gene.expr                <- tmp[tmp$Entrez_Gene_Id %in% pam50.gene.df$EntrezGene.ID,]
rownames(pam50.gene.expr)      <- pam50.gene.expr$Hugo_Symbol
annot.matrix                   <- cbind(pam50.gene.expr$Hugo_Symbol,pam50.gene.expr$Entrez_Gene_Id)
colnames(annot.matrix)         <- c('probe','EntrezGene.ID')
rownames(annot.matrix)         <- annot.matrix[,'probe']
pam50.gene.expr$Hugo_Symbol    <- NULL
pam50.gene.expr$Entrez_Gene_Id <- NULL
pam50.subtype.rs               <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr %>% t,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )

MBC.sample         <- colnames(pam50.gene.expr)
MBC.LumB.sample    <- MBC.sample[pam50.subtype.rs$subtype == 'LumB']
MBC.Basal.sample   <- MBC.sample[pam50.subtype.rs$subtype == 'Basal']
MBC.Her2.sample    <- MBC.sample[pam50.subtype.rs$subtype == 'Her2']
MBC.LumA.sample    <- MBC.sample[pam50.subtype.rs$subtype == 'LumA']
MBC.Normal.sample  <- MBC.sample[pam50.subtype.rs$subtype == 'Normal']


#####
load('client-side/output/organize.TCGA.and.MET500.breast.cancer.cnv.data.R.output/organize.TCGA.and.MET500.breast.cancer.cnv.data.RData')

tmp <- fread(input='client-side/Data/brca_mbcproject_wagle_2017/data_cna_hg19.seg') %>% as.data.frame()
require(CNTools)
MBC.seg                              <- CNSeg(tmp)
MBC.gene.cnv                         <- getRS(object = MBC.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "max")@rs
MBC.gene.cnv                         <- MBC.gene.cnv[,5:ncol(MBC.gene.cnv)]
rownames(MBC.gene.cnv)               <- MBC.gene.cnv$genename
MBC.gene.cnv$genename                <- NULL
MBC.gene.cnv.matrix                  <- as.matrix(MBC.gene.cnv)


draw.df <- MBC.sample.meta[MBC.Basal.sample,]
draw.df$cnv <- MBC.gene.cnv.matrix['NOTCH3',MBC.Basal.sample]
