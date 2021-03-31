require(GSVA)
require(dplyr)

gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id


G1.and.S.signature <- read.table("~/Project/BreastCancerMetaPotenial/client-side/Data/cell.cycle/G1.and.S.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
G2.and.M.signature <- read.table("~/Project/BreastCancerMetaPotenial/client-side/Data/cell.cycle/G2.and.M.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1


idx                            <- match(x=G1.and.S.signature,table = mapping.df$symbol)
G1.and.S.signature.ensemble.id <- mapping.df$ensembl_gene_id[idx]
G1.and.S.signature.ensemble.id <- G1.and.S.signature.ensemble.id[is.na(G1.and.S.signature.ensemble.id) == FALSE]
idx                            <- match(x=G2.and.M.signature,table = mapping.df$symbol)
G2.and.M.signature.ensemble.id <- mapping.df$ensembl_gene_id[idx]
G2.and.M.signature.ensemble.id <- G2.and.M.signature.ensemble.id[is.na(G2.and.M.signature.ensemble.id) == FALSE]


gm <- function(x) {
    x <- x[x > 0]  
    exp(mean(log(x)))
}

cell.cycle.score.list <- list()
perform.cell.cycle.analysis <- function(PRI.log2.tpm.matrix,MET.log2.tpm.matrix){
    v1                         <- apply(PRI.log2.tpm.matrix,1,function(x) median(x) > 1)
    v2                         <- apply(MET.log2.tpm.matrix,1,function(x) median(x) > 1)
    c.gene                     <- rownames(PRI.log2.tpm.matrix)[v1 & v2]
    PRI.score.matrix           <- gsva(expr = PRI.log2.tpm.matrix[c.gene,],gset.idx.list = list(G1.and.S = G1.and.S.signature.ensemble.id, G2.and.M = G2.and.M.signature.ensemble.id ),method='ssgsea',ssgsea.norm=FALSE)
    MET.score.matrix           <- gsva(expr = MET.log2.tpm.matrix[c.gene,],gset.idx.list = list(G1.and.S = G1.and.S.signature.ensemble.id, G2.and.M = G2.and.M.signature.ensemble.id ),method='ssgsea',ssgsea.norm=FALSE)
    G1.and.S.pri.score         <- PRI.score.matrix['G1.and.S',]
    G2.and.M.pri.score         <- PRI.score.matrix['G2.and.M',]
    G1.and.S.met.score         <- MET.score.matrix['G1.and.S',]
    G2.and.M.met.score         <- MET.score.matrix['G2.and.M',]
    
    list(G1.and.S.pri.score = G1.and.S.pri.score,G2.and.M.pri.score= G2.and.M.pri.score,G1.and.S.met.score = G1.and.S.met.score,G2.and.M.met.score = G2.and.M.met.score)  
}


################ Breast cancer ###################
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')


### Basal-like subtype #########
PRI.log2.tpm.matrix                   <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]
MET.log2.tpm.matrix                   <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.Basal.sample]
cell.cycle.score.list[['BRCA.Basal']] <- perform.cell.cycle.analysis(PRI.log2.tpm.matrix,MET.log2.tpm.matrix)

### LumB subtype #########
PRI.log2.tpm.matrix                  <- log2.tpm.matrix[,pure.PRI.breast.cancer.LumB.sample]
MET.log2.tpm.matrix                  <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.LumB.sample]
cell.cycle.score.list[['BRCA.LumB']] <- perform.cell.cycle.analysis(PRI.log2.tpm.matrix,MET.log2.tpm.matrix)


### Her2 subtype #########
PRI.log2.tpm.matrix                  <- log2.tpm.matrix[,pure.PRI.breast.cancer.Her2.sample]
MET.log2.tpm.matrix                  <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.Her2.sample]
cell.cycle.score.list[['BRCA.Her2']] <- perform.cell.cycle.analysis(PRI.log2.tpm.matrix,MET.log2.tpm.matrix)


############ Prostate cancer ########################3
load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')
load('server-side/RData//Prostate Adenocarcinoma.RData')
PRI.log2.tpm.matrix             <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]
MET.log2.tpm.matrix             <- MET500.log2.tpm.matrix[,pure.MET.prostate.cancer.sample]
cell.cycle.score.list[['PRAD']] <-  perform.cell.cycle.analysis(PRI.log2.tpm.matrix,MET.log2.tpm.matrix)


########### coleractal cancer ############
load('client-side/output/Select.pure.sample.colorectal.cancer.R.output/Select.pure.sample.colorectal.cancer.RData')
load('server-side/RData/COLORECTAL_SRP029880.RData')
PRI.log2.tpm.matrix             <- COLORECTAL_SRP029880_log2.tpm.matrix[,pure.PRI.colorectal.cancer.sample]
MET.log2.tpm.matrix             <- COLORECTAL_SRP029880_log2.tpm.matrix[,pure.MET.colorectal.cancer.sample]
cell.cycle.score.list[['COAD']] <- perform.cell.cycle.analysis(PRI.log2.tpm.matrix,MET.log2.tpm.matrix)


########## NET PAAD ############
load('client-side/output/Select.pure.sample.NET.pancreatic.cancer.R.output/Select.pure.sample.NET.pancreatic.cancer.RData')
load('server-side/RData/GEP.NET.RData')
PRI.log2.tpm.matrix             <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
MET.log2.tpm.matrix             <- GEP.NET.log2.tpm.matrix[,pure.MET.NET.pancreatic.cancer.sample]
cell.cycle.score.list[['PNET']] <- perform.cell.cycle.analysis(PRI.log2.tpm.matrix,MET.log2.tpm.matrix)


########## NET SI ############
load('client-side/output/Select.pure.sample.NET.si.cancer.R.output/Select.pure.sample.NET.si.cancer.RData')
load('server-side/RData/GEP.NET.RData')
PRI.log2.tpm.matrix              <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.si.cancer.sample]
MET.log2.tpm.matrix              <- GEP.NET.log2.tpm.matrix[,pure.MET.NET.si.cancer.sample]
cell.cycle.score.list[['SINET']] <- perform.cell.cycle.analysis(PRI.log2.tpm.matrix,MET.log2.tpm.matrix)


cell.cycle.score.df <- foreach(cancer = names(cell.cycle.score.list),.combine='rbind') %do% {
    o   <-   cell.cycle.score.list[[cancer]]
    df1 <- data.frame(phase = "G1/S",site='primary',    score =o$G1.and.S.pri.score)
    df2 <- data.frame(phase = "G2/M",site='primary',    score =o$G2.and.M.pri.score)
    df3 <- data.frame(phase = "G1/S",site='metastatic', score =o$G1.and.S.met.score)
    df4 <- data.frame(phase = "G2/M",site='metastatic', score =o$G2.and.M.met.score)
    
    tmp1        <- rbind(df1,df3)
    tmp1$score  <- scale(tmp1$score)
    tmp2        <- rbind(df2,df4)
    tmp2$score  <- scale(tmp2$score)
    
    df  <- rbind(tmp1,tmp2)
    df$cancer.type <- cancer
    df
}

#require(ggplot2)
#ggplot(cell.cycle.score.df[cell.cycle.score.df$phase == 'G1/S',], aes(x=cancer.type, y=score, fill=site)) + geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge(),size=3.0) + ylab('G1/S score')
#ggplot(cell.cycle.score.df[cell.cycle.score.df$phase == 'G2/M',], aes(x=cancer.type, y=score, fill=site)) + geom_boxplot(outlier.shape = NA) + geom_point(position=position_jitterdodge(),size=3.0) + ylab('G2/M score')

cell.cycle.score.df$cancer.type <- factor(cell.cycle.score.df$cancer.type,levels = c('BRCA.Basal','BRCA.LumB','BRCA.Her2','PRAD','COAD','PNET','SINET'))



phase.score.diff.df <- foreach(cancer = names(cell.cycle.score.list),.combine='rbind') %do% {
    o              <-   cell.cycle.score.list[[cancer]]
    df1            <- data.frame(site='primary',   score = o$G2.and.M.pri.score - o$G1.and.S.pri.score)
    df2            <- data.frame(site='metastatic',score = o$G2.and.M.met.score - o$G1.and.S.met.score)
    df             <- rbind(df1,df2)
    df$cancer.type <- cancer
    df$score       <- scale(df$score)
    df
}
#ggplot(phase.score.diff.df,aes(x=cancer.type,y=score,fill=site)) + geom_boxplot(outlier.shape = NA)  + geom_point(position=position_jitterdodge(),size=3.0)
phase.score.diff.df$cancer.type <- factor(phase.score.diff.df$cancer.type,levels = c('BRCA.Basal','BRCA.LumB','BRCA.Her2','PRAD','COAD','PNET','SINET'))



save(file='client-side/output/cell.cycle.R.output/cell.cycle.RData',list=c('cell.cycle.score.df','phase.score.diff.df'))




load('server-side/RData/PROSTATE_SRP253428.RData')
liver.met.id         <- PROSTATE_SRP253428_Metadata$Run[PROSTATE_SRP253428_Metadata$source_name == 'CRPC liver metastasis'] %>% as.character()
MET.log2.tpm.matrix  <- PROSTATE_SRP253428_log2.tpm.matrix[,liver.met.id]


load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')
load('server-side/RData//Prostate Adenocarcinoma.RData')
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]

score.rs                  <- perform.cell.cycle.analysis(PRI.log2.tpm.matrix,MET.log2.tpm.matrix)


G1.and.S.scroe.df <- rbind(
                     data.frame(score = score.rs$G1.and.S.pri.score,site = 'primary' ),
                     data.frame(score = score.rs$G1.and.S.met.score,site = 'metastatic')
                          )

G2.and.M.scroe.df <- rbind(
  data.frame(score = score.rs$G2.and.M.pri.score,site = 'primary' ),
  data.frame(score = score.rs$G2.and.M.met.score,site = 'metastatic' )
)

phase.diff.score.df <- data.frame(score = G2.and.M.scroe.df$score - G1.and.S.scroe.df$score,
                                  site = G1.and.S.scroe.df$site)

ggplot(G1.and.S.scroe.df, aes(x=site,y= scale(score), fill=site)) + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(),size=3.0,show.legend = FALSE) + 
  ylab('G1/S ssGSEA score') + ggplot.style  + scale_fill_manual(values = c('primary'='#EF8A62','metastatic'='#D1E5F0')) + ylim(-7,7)


ggplot(G2.and.M.scroe.df, aes(x=site,y= scale(score), fill=site)) + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(),size=3.0,show.legend = FALSE) + 
  ylab('G2/M ssGSEA score') + ggplot.style  + scale_fill_manual(values = c('primary'='#EF8A62','metastatic'='#D1E5F0')) + ylim(-7,7)


ggplot(phase.diff.score.df, aes(x=site,y= scale(score), fill=site)) + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(),size=3.0,show.legend = FALSE) + 
  ylab('G2/M - G1/S ssGSEA score') + ggplot.style  + scale_fill_manual(values = c('primary'='#EF8A62','metastatic'='#D1E5F0')) + ylim(-7,7)



save(file='client-side/output/cell.cycle.R.output/SRP253428.validation.RData',list=c('G1.and.S.scroe.df','G2.and.M.scroe.df','phase.diff.score.df'))






