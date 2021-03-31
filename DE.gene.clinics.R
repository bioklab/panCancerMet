library(readxl)
library(survival)
#library(ggsurvplot)
require(foreach)
require(GSVA)
require(dplyr)

## load TCGA clinical data ##
TCGA.CDR.data <- read_excel("client-side/Data/TCGA-CDR/TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR") %>% as.data.frame
get.patient.id <- function(x) { 
    tmp <- strsplit(x,split = '-') %>% unlist
    paste(tmp[1:3],collapse = '-')
}


## functions to perform the analysis ##
perform.suvival.analysis <- function(PRI.log2.tpm.matrix,target.gene.vec) {
  rs <- foreach(g = target.gene.vec,.combine='rbind')  %do%  {
        v                   <- PRI.log2.tpm.matrix[g,] 
        if( median(v) < log2(1 +1 )){
            return(data.frame(gene.name = g,hr=NA,p.val = NA))
        }
        high.group          <- sapply(names(v)[ v >= median(v)],get.patient.id)
        low.group           <- sapply(names(v)[ v <  median(v)],get.patient.id)
        high.group.df       <- TCGA.CDR.data[TCGA.CDR.data$bcr_patient_barcode %in% high.group,]
        low.group.df        <- TCGA.CDR.data[TCGA.CDR.data$bcr_patient_barcode %in% low.group,]
        high.group.df$group <- 'high'
        low.group.df$group  <- 'low'
        c.df                <- rbind(high.group.df,low.group.df)
        c.df$group          <- factor(c.df$group,levels = c('low','high'))
        sdf                 <- survdiff(Surv(PFI.time,PFI) ~  group, data = c.df)
        p.val               <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        hr                  <- log2((sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])) # according to https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
        data.frame(gene.name =g ,hr=hr,p.val = p.val)
  }
  rownames(rs) <- target.gene.vec
  rs[complete.cases(rs),]
}

perform.stage.analysis <- function(PRI.log2.tpm.matrix,up.gene,dn.gene) {
    v                   <- apply(PRI.log2.tpm.matrix,1,function(x) median(x) > 1)
    PRI.log2.tpm.matrix <- PRI.log2.tpm.matrix[v,]
    ssgsea.score.matrix <- gsva(expr = PRI.log2.tpm.matrix,
                                method='ssgsea',
                                ssgsea.norm=FALSE,
                                gset.idx.list=list(up.gene = up.gene,dn.gene = dn.gene)
    )
    patient.id <- sapply(PRI.log2.tpm.matrix %>% colnames,get.patient.id)
    df <- data.frame(up.score = ssgsea.score.matrix['up.gene',],
                     dn.score = ssgsea.score.matrix['dn.gene',],
                     stage    = TCGA.CDR.data$ajcc_pathologic_tumor_stage [match(patient.id,TCGA.CDR.data$bcr_patient_barcode)]     
                     )
    df <- df[complete.cases(df),]
    mapping.vec <- c('Stage I' = 'Stage I', 'Stage IA' = 'Stage I','Stage IA' = 'Stage IB',
                     'Stage II' = 'Stage II', 'Stage IIA' = 'Stage II','Stage IIB' = 'Stage II','Stage IIC' = 'Stage II',
                     'Stage III' = 'Stage II', 'Stage IIIA' = 'Stage III','Stage IIIB' = 'Stage III','Stage IIIC' = 'Stage III'
                     )

    df$stage <- mapping.vec[df$stage %>% as.character()]
    df[complete.cases(df),]
}

perform.stage.analysis.for.single.gene <- function(PRI.log2.tpm.matrix,up.gene,dn.gene) {
    expr.matrix <- PRI.log2.tpm.matrix[c(up.gene,dn.gene),]
    patient.id  <- sapply(PRI.log2.tpm.matrix %>% colnames,get.patient.id)
    stage       <- TCGA.CDR.data$ajcc_pathologic_tumor_stage [match(patient.id,TCGA.CDR.data$bcr_patient_barcode)]  
    mapping.vec <- c('Stage I' = 'Stage I', 'Stage IA' = 'Stage I','Stage IA' = 'Stage IB',
                     'Stage II' = 'Stage II', 'Stage IIA' = 'Stage II','Stage IIB' = 'Stage II','Stage IIC' = 'Stage II',
                     'Stage III' = 'Stage II', 'Stage IIIA' = 'Stage III','Stage IIIB' = 'Stage III','Stage IIIC' = 'Stage III'
    )
    stage    <- mapping.vec[stage]
    s1.flag <- stage == 'Stage I'
    s2.flag <- stage == 'Stage II'
    s3.flag <- stage == 'Stage III'
    
    
    up.p.value.df <- foreach(g = up.gene,.combine='rbind') %do% {
       s1.expr <- PRI.log2.tpm.matrix[g,s1.flag]  
       s2.expr <- PRI.log2.tpm.matrix[g,s2.flag]  
       s3.expr <- PRI.log2.tpm.matrix[g,s3.flag]  
       data.frame(II.vs.I.p.value   = wilcox.test(s1.expr,s2.expr)$p.value,
                  III.vs.II.p.value = wilcox.test(s2.expr,s3.expr)$p.value)
    }
    up.p.value.df$gene <- up.gene
    up.p.value.df$type <- 'up'
    
    dn.p.value.df <- foreach(g = dn.gene,.combine='rbind') %do% {
      s1.expr <- PRI.log2.tpm.matrix[g,s1.flag]  
      s2.expr <- PRI.log2.tpm.matrix[g,s2.flag]  
      s3.expr <- PRI.log2.tpm.matrix[g,s3.flag]  
      data.frame(II.vs.I.p.value = wilcox.test(s1.expr,s2.expr)$p.value,
                 III.vs.II.p.value = wilcox.test(s2.expr,s3.expr)$p.value)
    }
    dn.p.value.df$gene <- dn.gene
    dn.p.value.df$type <- 'dn'
    rbind(up.p.value.df,dn.p.value.df)
}




## let us perform analysis for BRCA ##
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')

## BRCA.Basal ##
PRI.log2.tpm.matrix    <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]
DE.rs                  <- BRCA.Basal.DE.rs
up.df                  <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$up.gene)
dn.df                  <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
up.df$gene.type        <- 'up'
dn.df$gene.type        <- 'dn'

BRCA.Basal.survival.df             <- rbind(up.df,dn.df)
BRCA.Basal.survival.df$cancer.type <- 'BRCA.Basal'

BRCA.Basal.stage.df                <- perform.stage.analysis(PRI.log2.tpm.matrix,DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
BRCA.Basal.stage.df$cancer.type    <- 'BRCA.Basal'
BRCA.Basal.stage.df$up.score       <- scale(BRCA.Basal.stage.df$up.score)
BRCA.Basal.stage.df$dn.score       <- scale(BRCA.Basal.stage.df$dn.score)

BRCA.Basal.stage.single.gene.df                <- perform.stage.analysis.for.single.gene(PRI.log2.tpm.matrix,DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
BRCA.Basal.stage.single.gene.df$cancer.type    <- 'BRCA.Basal'

fdr1 <- p.adjust(BRCA.Basal.stage.single.gene.df$II.vs.I.p.value,method='fdr')
fdr2 <- p.adjust(BRCA.Basal.stage.single.gene.df$III.vs.II.p.value,method='fdr')
sum(fdr1 < 0.05)
sum(fdr2 < 0.05)


## BRCA.LumB ##
PRI.log2.tpm.matrix    <- log2.tpm.matrix[,pure.PRI.breast.cancer.LumB.sample]
DE.rs                  <- BRCA.LumB.DE.rs
up.df                  <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$up.gene)
dn.df                  <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
up.df$gene.type        <- 'up'
dn.df$gene.type        <- 'dn'

BRCA.LumB.survival.df             <- rbind(up.df,dn.df)
BRCA.LumB.survival.df$cancer.type <- 'BRCA.LumB'

BRCA.LumB.stage.df                <- perform.stage.analysis(PRI.log2.tpm.matrix,DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
BRCA.LumB.stage.df$cancer.type    <- 'BRCA.LumB'
BRCA.LumB.stage.df$up.score       <- scale(BRCA.LumB.stage.df$up.score)
BRCA.LumB.stage.df$dn.score       <- scale(BRCA.LumB.stage.df$dn.score)

BRCA.LumB.stage.single.gene.df                <- perform.stage.analysis.for.single.gene(PRI.log2.tpm.matrix,DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
BRCA.LumB.stage.single.gene.df$cancer.type    <- 'BRCA.LumB'

fdr1 <- p.adjust(BRCA.LumB.stage.single.gene.df$II.vs.I.p.value,method='fdr')
fdr2 <- p.adjust(BRCA.LumB.stage.single.gene.df$III.vs.II.p.value,method='fdr')
sum(fdr1 < 0.05)
sum(fdr2 < 0.05)

## BRCA.Her2 ##
PRI.log2.tpm.matrix    <- log2.tpm.matrix[,pure.PRI.breast.cancer.Her2.sample]
DE.rs                  <- BRCA.Her2.DE.rs
up.df                  <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$up.gene)
dn.df                  <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
up.df$gene.type        <- 'up'
dn.df$gene.type        <- 'dn'

BRCA.Her2.survival.df             <- rbind(up.df,dn.df)
BRCA.Her2.survival.df$cancer.type <- 'BRCA.Her2'

BRCA.Her2.stage.df                <- perform.stage.analysis(PRI.log2.tpm.matrix,DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
BRCA.Her2.stage.df$cancer.type    <- 'BRCA.Her2'
BRCA.Her2.stage.df$up.score       <- scale(BRCA.Her2.stage.df$up.score)
BRCA.Her2.stage.df$dn.score       <- scale(BRCA.Her2.stage.df$dn.score)

BRCA.Her2.stage.single.gene.df                <- perform.stage.analysis.for.single.gene(PRI.log2.tpm.matrix,DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
BRCA.Her2.stage.single.gene.df$cancer.type    <- 'BRCA.Her2'

fdr1 <- p.adjust(BRCA.Her2.stage.single.gene.df$II.vs.I.p.value,method='fdr')
fdr2 <- p.adjust(BRCA.Her2.stage.single.gene.df$III.vs.II.p.value,method='fdr')
sum(fdr1 < 0.05)
sum(fdr2 < 0.05)

## let us perform analysis for PRAD ##
load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')
load('server-side/RData//Prostate Adenocarcinoma.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')

PRI.log2.tpm.matrix          <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]
DE.rs                        <- PRAD.DE.rs
up.df                        <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$up.gene)
dn.df                        <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)

up.df$gene.type              <- 'up'
dn.df$gene.type              <- 'dn'
PRAD.survival.df             <- rbind(up.df,dn.df)
PRAD.survival.df$cancer.type <- 'PRAD'

#PRAD, stage information is NOT available


## let us perform analysis for COAD ##
load('server-side/RData//Colon Adenocarcinoma.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')

PRI.log2.tpm.matrix        <- log2.tpm.matrix
DE.rs                      <- COAD.DE.rs
up.df                      <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$up.gene)
dn.df                      <- perform.suvival.analysis(PRI.log2.tpm.matrix,   DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
up.df$gene.type            <- 'up'
dn.df$gene.type            <- 'dn'

COAD.survival.df             <- rbind(up.df,dn.df)
COAD.survival.df$cancer.type <- 'COAD'

COAD.stage.df                <- perform.stage.analysis(PRI.log2.tpm.matrix,DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
COAD.stage.df$cancer.type    <- 'COAD'
COAD.stage.df$up.score       <- scale(COAD.stage.df$up.score)
COAD.stage.df$dn.score       <- scale(COAD.stage.df$dn.score)

COAD.stage.single.gene.df                <- perform.stage.analysis.for.single.gene(PRI.log2.tpm.matrix,DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene)
COAD.stage.single.gene.df$cancer.type    <- 'COAD'

fdr1 <- p.adjust(COAD.stage.single.gene.df$II.vs.I.p.value,method='fdr')
fdr2 <- p.adjust(COAD.stage.single.gene.df$III.vs.II.p.value,method='fdr')
sum(fdr1 < 0.05)
sum(fdr2 < 0.05)


## ok, save results ##
survival.df <- rbind(BRCA.Basal.survival.df,BRCA.LumB.survival.df,BRCA.Her2.survival.df,PRAD.survival.df,COAD.survival.df)
stage.df    <- rbind(BRCA.Basal.stage.df,BRCA.LumB.stage.df,BRCA.Her2.stage.df,COAD.stage.df)

save(file='client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData',list=c('survival.df','stage.df'))


# cc <- rbind(BRCA.Basal.survival.df,BRCA.LumB.survival.df,BRCA.Her2.survival.df,PRAD.survival.df,COAD.survival.df)
# 
# up.g <- cc[cc$hr >0 &cc$p.val < 0.05 & cc$gene.type == 'up',]
# 
# dn.g <- cc[cc$hr <0 &cc$p.val < 0.05 & cc$gene.type == 'dn',]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# load('client-side/output/Select.pure.sample.colorectal.cancer.R.output/Select.pure.sample.colorectal.cancer.RData')
# load('server-side/RData/COLORECTAL_SRP029880.RData')
# PRI.log2.read.count.matrix <- COLORECTAL_SRP029880_log2.read.count.matrix[,pure.PRI.colorectal.cancer.sample]
# PRI.log2.tpm.matrix        <- COLORECTAL_SRP029880_log2.tpm.matrix[,pure.PRI.colorectal.cancer.sample]
# #PRRX1.DE.rs.list[['COAD']]       <- perform.DE.between.PRRX1.high.and.low()
# 
# 
# MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.MET.colorectal.cancer.sample]
# df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
# df$disease <- 'COAD'
# PRRX1.expr.rs.list[['COAD']] <- df
# 
# 
# 
# 
# PRI.log2.tpm.matrix        <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]
# #PRRX1.DE.rs.list[['PRAD']]       <- perform.DE.between.PRRX1.high.and.low()
# 
# 
# MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.MET.prostate.cancer.sample]
# df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
# df$disease <- 'PRAD'
# PRRX1.expr.rs.list[['PRAD']] <- df
# 
# 
# load('client-side/output/Select.pure.sample.NET.pancreatic.cancer.R.output/Select.pure.sample.NET.pancreatic.cancer.RData')
# load('server-side/RData/GEP.NET.RData')
# PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
# PRI.log2.tpm.matrix       <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
# #PRRX1.DE.rs.list[['NET.PAAD']]   <- perform.DE.between.PRRX1.high.and.low()
# 
# 
# MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
# df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
# df$disease <- 'NET.PAAD'
# PRRX1.expr.rs.list[['NET.PAAD']] <- df
# 
# 
# 
# load('client-side/output/Select.pure.sample.NET.si.cancer.R.output/Select.pure.sample.NET.si.cancer.RData')
# load('server-side/RData/GEP.NET.RData')
# PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.si.cancer.sample]
# PRI.log2.tpm.matrix        <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.si.cancer.sample]
# #PRRX1.DE.rs.list[['NET.SI']]     <- perform.DE.between.PRRX1.high.and.low()
# 
# MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.MET.NET.si.cancer.sample]
# df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
# df$disease <- 'NET.SI'
# PRRX1.expr.rs.list[['NET.SI']] <- df
# 
# 
# 
# 
# 
# patient.id <- sapply(pure.PRI.prostate.cancer.sample,get.patient.id)
# df         <- TCGA.CDR.data[TCGA.CDR.data$bcr_patient_barcode %in% patient.id,]
# df         <- df[,c('PFI.time','PFI')]
# 
# fit <- survfit(Surv(PFI.time,PFI) ~ 1, data = df)





