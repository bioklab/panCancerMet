require(data.table)
require(reshape2)

get.gene.name <- function(x) {
    strsplit(x=x,split = "\\.\\.")[[1]][1]  
}

############### organize CRISPR screnn results ########################
Achilles_gene_effect                  <- read.csv("client-side/Data/DepMap/Achilles_gene_effect_unscaled.csv", stringsAsFactors=FALSE)
rownames(Achilles_gene_effect)        <- Achilles_gene_effect$X
Achilles_gene_effect$X                <- NULL
Achilles_gene_effect.matrix           <- as.matrix(Achilles_gene_effect)
tmp                                   <- colnames(Achilles_gene_effect.matrix)
name.vec                              <- sapply(tmp,get.gene.name)
names(name.vec)                       <- NULL
colnames(Achilles_gene_effect.matrix) <- name.vec
Achilles.gene.effect.matrix           <- Achilles_gene_effect.matrix

# sample_info <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/sample_info.csv", stringsAsFactors=FALSE)
# rownames(sample_info) <- sample_info$DepMap_ID
# rownames(Achilles_gene_effect.matrix) <- sample_info[rownames(Achilles_gene_effect.matrix),'CCLE.Name']



############### organize RNAi screnn results ########################
RNAi.df <- fread(input='client-side/Data/DepMap/D2_combined_gene_dep_scores.csv',header = TRUE) %>% as.data.frame()
get.gene.name.v2 <- function(x) {
  strsplit(x=x,split = " ")[[1]][1]  
}
tmp       <- sapply(RNAi.df$V1, get.gene.name.v2)
names(tmp) <- NULL

RNAi.df$V1 <- NULL
RNAi.screen.matrix <- as.matrix(RNAi.df)
rownames(RNAi.screen.matrix) <- tmp

############# Cell line information ############
sample.info           <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/sample_info.csv", stringsAsFactors=FALSE)
rownames(sample.info) <- sample.info$DepMap_ID
CRISPR.screen.cell.line.info  <- sample.info



############### organize CCLE cnv results ########################

# CCLE.gene.cnv           <- read.csv("client-side/Data/DepMap/CCLE_gene_cn.csv", stringsAsFactors=FALSE)
# rownames(CCLE.gene.cnv) <- CCLE.gene.cnv$X
# CCLE.gene.cnv$X         <- NULL
# CCLE.gene.cnv.matrix    <- as.matrix(CCLE.gene.cnv)
# tmp                            <- colnames(CCLE.gene.cnv.matrix)
# name.vec                       <- sapply(tmp,get.gene.name)
# names(name.vec)                <- NULL
# colnames(CCLE.gene.cnv.matrix) <- name.vec

# sample_info <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/sample_info.csv", stringsAsFactors=FALSE)
# rownames(sample_info) <- sample_info$DepMap_ID
# rownames(Achilles_gene_effect.matrix) <- sample_info[rownames(Achilles_gene_effect.matrix),'CCLE.Name']


############### organize drug sensitivity results ########################
# require(data.table)
# primary.drug.screen              <- fread(input = "client-side/Data/DepMap/primary-screen-replicate-collapsed-logfold-change.csv", stringsAsFactors=FALSE) %>% as.data.frame
# rownames(primary.drug.screen)    <- primary.drug.screen$V1
# primary.drug.screen$V1           <- NULL
# primary.drug.screen.matrix       <- as.matrix(primary.drug.screen)
# 
# secondary.drug.screen            <- fread(input = "client-side/Data/DepMap/secondary-screen-replicate-collapsed-logfold-change.csv", stringsAsFactors=FALSE) %>% as.data.frame
# rownames(secondary.drug.screen)  <- secondary.drug.screen$V1
# secondary.drug.screen$V1         <- NULL
# secondary.drug.screen.matrix     <- as.matrix(secondary.drug.screen)
# 
# primary.drug.screen.meta               <- fread(input='client-side/Data/DepMap/primary-screen-replicate-collapsed-treatment-info.csv',   stringsAsFactors=FALSE) %>% as.data.frame
# secondary.drug.screen.meta             <- fread(input='client-side/Data/DepMap/secondary-screen-replicate-collapsed-treatment-info.csv', stringsAsFactors=FALSE) %>% as.data.frame
# rownames(primary.drug.screen.meta)     <- primary.drug.screen.meta$column_name
# primary.drug.screen.meta$column_name   <- NULL
# rownames(secondary.drug.screen.meta)   <- secondary.drug.screen.meta$column_name
# secondary.drug.screen.meta$column_name <- NULL

################## CCLE expression ##############
# require(data.table)
# CCLE_expression <- read.csv("client-side/Data/DepMap/CCLE_expression.csv", stringsAsFactors=FALSE) 
# rownames(CCLE_expression) <- CCLE_expression$X
# CCLE_expression$X        <- NULL
# CCLE_expression.matrix    <- as.matrix(CCLE_expression)
# tmp                  <- colnames(CCLE_expression.matrix)
# get.gene.name <- function(x) {
#   strsplit(x=x,split = "\\.\\.")[[1]][1]  
# }
# name.vec                       <- sapply(tmp,get.gene.name)
# names(name.vec)                <- NULL
# colnames(CCLE_expression.matrix) <- name.vec
# CCLE.expression.matrix <- CCLE_expression.matrix





############## Organize GDSC AUC dataset ###########
GDSC1_fitted_dose_response_15Oct19   <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/GDSC/GDSC1_fitted_dose_response_15Oct19.csv", stringsAsFactors=FALSE)
GDSC2_fitted_dose_response_15Oct19   <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/GDSC/GDSC2_fitted_dose_response_15Oct19.csv", stringsAsFactors=FALSE)
df                                   <- rbind(GDSC1_fitted_dose_response_15Oct19,GDSC2_fitted_dose_response_15Oct19)
df                                   <- df[,c('CELL_LINE_NAME','DRUG_NAME','AUC')]

GDSC.matrix               <- dcast(df,formula =  DRUG_NAME ~ CELL_LINE_NAME,value.var = 'AUC',fun.aggregate = mean)
rownames(GDSC.matrix)     <- GDSC.matrix$DRUG_NAME
GDSC.matrix$DRUG_NAME     <- NULL
GDSC.AUC.matrix           <- as.matrix(GDSC.matrix)



############## Organize GDSC IC50 dataset ###########
GDSC1_fitted_dose_response_15Oct19   <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/GDSC/GDSC1_fitted_dose_response_15Oct19.csv", stringsAsFactors=FALSE)
GDSC2_fitted_dose_response_15Oct19   <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/GDSC/GDSC2_fitted_dose_response_15Oct19.csv", stringsAsFactors=FALSE)
df                                   <- rbind(GDSC1_fitted_dose_response_15Oct19,GDSC2_fitted_dose_response_15Oct19)
df                                   <- df[,c('CELL_LINE_NAME','DRUG_NAME','LN_IC50')]

GDSC.matrix               <- dcast(df,formula =  DRUG_NAME ~ CELL_LINE_NAME,value.var = 'LN_IC50',fun.aggregate = mean)
rownames(GDSC.matrix)     <- GDSC.matrix$DRUG_NAME
GDSC.matrix$DRUG_NAME     <- NULL
GDSC.IC50.matrix           <- as.matrix(GDSC.matrix)



############## Organize CTRP AUC dataset ###########

AUC <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/CTRP/AUC.csv", stringsAsFactors=FALSE)
cell.line.annotation <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/CTRP/cell.line.annotation.csv", stringsAsFactors=FALSE)
cp.annotation <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/Data/CTRP/cp.annotation.csv", stringsAsFactors=FALSE)

cell.line.annotation <- cell.line.annotation[,c('ccl_name','master_ccl_id')]
cp.annotation        <- cp.annotation[,c('cpd_name','master_cpd_id')]


AUC <- AUC[,2:4]

tmp1 <- merge(x=AUC,y=cell.line.annotation,by='master_ccl_id',all.x=TRUE)
tmp2 <- merge(x=tmp1,y=cp.annotation,by = 'master_cpd_id',all.x=TRUE)
AUC.df <- tmp2[,c('ccl_name','cpd_name','area_under_curve')]


CTRP.AUC.matrix       <- dcast(AUC.df,formula =  cpd_name ~ ccl_name,value.var = 'area_under_curve',fun.aggregate = mean)
rownames(CTRP.AUC.matrix) <- CTRP.AUC.matrix$cpd_name
CTRP.AUC.matrix$cpd_name <- NULL
CTRP.AUC.matrix <- as.matrix(CTRP.AUC.matrix)

# df                      <- rbind(GDSC1_fitted_dose_response_15Oct19,GDSC2_fitted_dose_response_15Oct19)
# breast.cancer.cell.line <- df$CELL_LINE_NAME[df$TCGA_DESC == 'BRCA'] %>% unique
# breast.GDSC.matrix      <- GDSC.matrix[,breast.cancer.cell.line]
# 
# rename.cell.line <- function(x) {
#     x <- gsub(pattern = '-',x = x,replacement = '')  
#     toupper(x)
# }
# tmp                          <- sapply(colnames(breast.GDSC.matrix),rename.cell.line)
# names(tmp)                   <- NULL
# tmp                          <- paste(tmp,'_BREAST',sep='')
# colnames(breast.GDSC.matrix) <- tmp
# breast.GDSC.matrix           <- t(breast.GDSC.matrix)

############## Save data ###########
#save(list = c('Achilles.gene.effect.matrix','CCLE.expression.matrix','primary.drug.screen.matrix','primary.drug.screen.meta','secondary.drug.screen.matrix','secondary.drug.screen.meta','sample.info','CCLE.gene.cnv.matrix','breast.GDSC.matrix'),file = 'client-side/output/organize.Depmap.data.R.output/organize.Depmap.data.RData')
save(list = c('Achilles.gene.effect.matrix','RNAi.screen.matrix','CRISPR.screen.cell.line.info','GDSC.AUC.matrix','GDSC.IC50.matrix','CTRP.AUC.matrix'),file = 'client-side/output/organize.Depmap.data.R.output/organize.Depmap.data.RData')









###################
# c.cell.line    <- intersect(rownames(CCLE.gene.cnv.matrix),rownames(Achilles_gene_effect.matrix))
# tmp            <- cor(CCLE.gene.cnv.matrix[c.cell.line,'NOTCH3'],Achilles_gene_effect.matrix[c.cell.line,],method='spearman')
# cor.vec        <- c(tmp)
# names(cor.vec) <- colnames(tmp)
# cor.vec <- sort(cor.vec,decreasing = TRUE)
# 
# 
# 
# c.cell.line    <- intersect(rownames(CCLE.gene.cnv.matrix),rownames(CCLE_expression.matrix))
# tmp            <- cor(CCLE.gene.cnv.matrix[c.cell.line,'NOTCH3'],CCLE_expression.matrix[c.cell.line,],method='spearman')
# cor.vec        <- c(tmp)
# names(cor.vec) <- colnames(tmp)
# cor.vec <- sort(cor.vec,decreasing = TRUE)
# 
# 
# c.cell.line    <- intersect(rownames(CCLE_expression.matrix),rownames(CCLE_expression.matrix))
# tmp            <- cor(CCLE_expression.matrix[c.cell.line,'NOTCH3'],CCLE_expression.matrix[c.cell.line,],method='spearman')
# cor.vec        <- c(tmp)
# names(cor.vec) <- colnames(tmp)
# cor.vec <- sort(cor.vec,decreasing = TRUE)


