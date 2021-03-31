frs <- foreach(s = ME500.brca.sample) %do% {
  met.expr        <- MET500.log2.fpkm.matrix[liver.specific.gene,s]
  df              <- data.frame(liver.expr = liver.expr,met.expr=met.expr)
  flag            <- df$liver.expr >0 & df$met.expr > 0
  df              <- df[flag,]
  df$met.expr     <- log2(2^df$met.expr - 1)
  df$liver.expr   <- log2(2^df$liver.expr - 1)
  intercept       <- median(df$met.expr - df$liver.expr )
  df$residual     <- df$met.expr - (df$liver.expr + intercept)
  df$color <- 'OTHER'
  df[df %>% rownames == DHODH,'color'] <- 'DHODH'
  ggplot(df,aes(x=liver.expr,y=met.expr,color=color)) + geom_point() + geom_abline(slope=1,intercept = intercept) 
  
}



require(plyr)
require(dplyr)
require(genefu)
require(Rtsne)
require(ggplot2)
require(dplyr)
source('client-side/code/util.R')
load('server-side/RData/Breast Invasive Carcinoma.RData')


load('server-side/RData/BRACA_SRP043470.RData')

#### Perform pam50 subtyping for TCGA polyA samples###########
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#') # well, I borrow some information from the metastatic breast cancer evaluation project 
pam50.gene                  <- pam50.gene.df$ensemble.gene.id %>% as.character
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)

m.sample         <- BRACA_SRP043470_Metadata$Run[BRACA_SRP043470_Metadata$source_name=='Liver Metastasis Tumour'] %>% as.character()
pam50.gene.expr             <- BRACA_SRP043470_log2.fpkm.matrix[pam50.gene.df$probe %>% as.character,m.sample] %>% t 

m.sample <- c('SRR5357752','SRR5357760','SRR5357767') #From project SRP102119
pam50.gene.expr             <- LIVERMET.log2.fpkm.matrix[pam50.gene.df$probe %>% as.character,m.sample] %>% t 

m.sample <- c('SRR6298032','SRR6298035','SRR6298037','SRR6298040','SRR6298041') # TNBC
pam50.gene.expr             <- LIVERMET.log2.fpkm.matrix[pam50.gene.df$probe %>% as.character,m.sample] %>% t 



annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )
