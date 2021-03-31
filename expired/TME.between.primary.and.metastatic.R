#library(xCell)
library(plyr)
library(dplyr)
library(foreach)
require(ggplot2)
source('client-side/code/Manuscript/ggplot.style.R')


####################### prepare cancer data for analysis ###########################
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('client-side/output/immune.cell.marker.gene.R.output/immune.cell.marker.gene.RData')
load('client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')
load('client-side/output/estimate.hepatocytes.abundance.R.output/estimate.hepatocytes.abundance.RData')

TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
MET500.liver.sample                       <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
cancer.data                               <- cbind(MET500.log2.fpkm.matrix,TCGA.breast.cancer.log2.fpkm.matrix[rownames(MET500.log2.fpkm.matrix),])





############### function to compute immune cell score with signatures ########
compute.immune.score  <- function(.immune.cell.signature){
    flag                   <- .immune.cell.signature$ensemble.id %in% rownames(cancer.data)
    .immune.cell.signature <- .immune.cell.signature[flag,]
    cell.type.vec                         <- unique(.immune.cell.signature$cell.type) %>% as.character()
    immune.score.matrix <- foreach(cell.type = cell.type.vec,.combine='rbind') %do% {
        g <- .immune.cell.signature$ensemble.id[.immune.cell.signature$cell.type == cell.type] %>% as.character  
        if(length(g) > 1){
            x <- apply(cancer.data[g,],2,function(x) 2^mean(x))
            c(x)
         }
        else{
            2^cancer.data[g,] 
         }
    }
    if(length(cell.type.vec) == 1) {
        m           <- matrix(c(immune.score.matrix),nrow=1)   
        rownames(m) <- cell.type.vec
        colnames(m) <- colnames(cancer.data)
        immune.score.matrix <- m
    }else{
        rownames(immune.score.matrix) <- cell.type.vec
    }
    immune.score.matrix
}

######################## function for immune score analysis ##############################
compare.TME <- function(.case.sample,.control.sample,.immune.score.matrix,.purity.scale=FALSE,fit.method='loess'){
    purity.case      <- tumor.purity.based.on.cell.line.vec[.case.sample]
    purity.control   <- tumor.purity.based.on.cell.line.vec[.control.sample]
    if(.purity.scale == TRUE){
        purity.case           <- scale(purity.case,center = TRUE,scale = TRUE) %>% c
        purity.control        <- scale(purity.control,center=TRUE,scale=TRUE) %>% c
        names(purity.case)    <- .case.sample
        names(purity.control) <- .control.sample
        
    }
    x                <- c(purity.case,purity.control)
    cell.type.vec    <- rownames(.immune.score.matrix)
    
    rs <- foreach(cell.type = cell.type.vec,.combine='rbind') %do% {
        df                <- data.frame(x=x,y=.immune.score.matrix[cell.type,names(x)])
        df$y              <- log2(df$y)
        q.vec             <- quantile(df$y)
        cut.off           <- q.vec[4] + 1.5 * IQR(df$y)
        df                <- df[df$y < cut.off,]
        cut.off           <- q.vec[2] - 1.5 * IQR(df$y)
        df                <- df[df$y > cut.off,]
        if(fit.method == 'lm'){
            lm.fit.rs        <-  lm(data = df,formula= y ~ x) 
            adjusted.score   <- lm.fit.rs$residuals
        }else{
            loess.fit.rs     <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
            adjusted.score   <- loess.fit.rs$residuals
        }
        ..case.sample     <- intersect(.case.sample,rownames(df))
        ..control.sample  <- intersect(.control.sample,rownames(df))
        
       
        p.value           <- wilcox.test(adjusted.score[..case.sample],adjusted.score[..control.sample])$p.value
        delta             <- median(adjusted.score[..case.sample]) -  median(adjusted.score[..control.sample])
        data.frame(effect.size=delta,p.value=p.value)
    }
    rownames(rs) <- cell.type.vec  
    rs$fdr       <- p.adjust(rs$p.value,method='fdr')
    
    
    
    data.list <- foreach(cell.type = cell.type.vec) %do% {
        df                <- data.frame(x=x,y=.immune.score.matrix[cell.type,names(x)])
        df$y              <- log2(df$y)
        q.vec             <- quantile(df$y)
        cut.off           <- q.vec[4] + 1.5 * IQR(df$y)
        df                <- df[df$y < cut.off,]
        cut.off           <- q.vec[2] - 1.5 * IQR(df$y)
        df                <- df[df$y > cut.off,]
        if(fit.method == 'lm'){
            lm.fit.rs <-  lm(data = df,formula= y ~ x) 
            adjusted.score     <- lm.fit.rs$residuals
        }else{
            loess.fit.rs      <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
            adjusted.score     <- loess.fit.rs$residuals
        }
        ..case.sample     <- intersect(.case.sample,rownames(df))
        ..control.sample  <- intersect(.control.sample,rownames(df))
        
        case.df           <- data.frame(immune.score          = .immune.score.matrix[cell.type,..case.sample] %>% log2,
                                        adjusted.immune.score = adjusted.score[..case.sample],
                                        purity                = tumor.purity.based.on.cell.line.vec[..case.sample]  ,
                                        class                 = 'MET500')
        
        control.df           <- data.frame(immune.score       = .immune.score.matrix[cell.type,..control.sample] %>% log2,
                                        adjusted.immune.score = adjusted.score[..control.sample],
                                        purity                = tumor.purity.based.on.cell.line.vec[..control.sample] ,
                                        class                 = 'TCGA')
        rbind(case.df,control.df)
        
        
        
        
        # p.value           <- wilcox.test(adjusted.expr[..case.sample],adjusted.expr[..control.sample])$p.value
        # delta             <- median(adjusted.expr[..case.sample]) -  median(adjusted.expr[..control.sample])
        # data.frame(effect.size=delta,p.value=p.value)
    }
    names(data.list) <- cell.type.vec
    #rownames(rs) <- cell.type.vec  
    #rs$fdr       <- p.adjust(rs$p.value,method='fdr')
    
    # figure.list <- foreach(cell.type = cell.type.vec) %do% {
    #     df                <- data.frame(x=x,y=.immune.score.matrix[cell.type,names(x)])
    #     df$data.source    <- ifelse(rownames(df) %in% .control.sample,'control','case')
    #     p                 <- ggplot(df,aes(x=x,y=y)) + geom_point(aes(color=factor(data.source)),size=2) + geom_smooth(method=fit.method,se=FALSE) + ggplot.style +  theme(legend.position = "none") + xlab('purity') + ylab('immune.score') + ggtitle(cell.type)
    #     p
    # }
    # names(figure.list) <- cell.type.vec
    # 
    # data.list <- foreach(cell.type = cell.type.vec) %do% {
    #     df                <- data.frame(x=x,y=.immune.score.matrix[cell.type,names(x)])
    #     df$data.source    <- ifelse(rownames(df) %in% .control.sample,'control','case')
    #     df
    # }
    # names(data.list) <- cell.type.vec
    
    #list(rs=rs,figures=figure.list,data=data.list)
    list(rs=rs,data=data.list)
    
}

######### cell types considered in this study ##########
cell.type.vec <- immune.cell.signature$cell.type %>% unique


#################  Liver metastasis vs primary tumor, subtype-specific TME comaprision  ########################

immune.score.matrix  <- compute.immune.score(immune.cell.signature)

MET500.sample        <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Basal.sample
M.vs.P.Basal.rs      <- compare.TME(.case.sample = MET500.sample,.control.sample = TCGA.sample,.immune.score.matrix = immune.score.matrix)
cell.type.vec[M.vs.P.Basal.rs$rs$fdr < 0.01]
cell.type.vec[M.vs.P.Basal.rs$rs$fdr < 0.05]


MET500.sample        <- MET500.breast.cancer.polyA.LumB.sample
MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample          <- pure.TCGA.breast.cancer.polyA.LumB.sample
M.vs.P.LumB.rs       <- compare.TME(.case.sample = MET500.sample,.control.sample = TCGA.sample,.immune.score.matrix = immune.score.matrix)
cell.type.vec[M.vs.P.LumB.rs$rs$fdr < 0.01]
cell.type.vec[M.vs.P.LumB.rs$rs$fdr < 0.05]


MET500.sample        <- MET500.breast.cancer.polyA.Her2.sample
MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Her2.sample
M.vs.P.Her2.rs       <- compare.TME(.case.sample = MET500.sample,.control.sample = TCGA.sample,.immune.score.matrix = immune.score.matrix)
cell.type.vec[M.vs.P.Her2.rs$rs$fdr < 0.01]
cell.type.vec[M.vs.P.Her2.rs$rs$fdr < 0.05]


# MET500.sample        <- MET500.breast.cancer.polyA.LumA.sample
# MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
# MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
# TCGA.sample          <- pure.TCGA.breast.cancer.polyA.LumA.sample
# M.vs.P.LumA.rs       <- compare.TME(.case.sample = MET500.sample,.control.sample = TCGA.sample,.immune.score.matrix = immune.score.matrix)
# cell.type.vec[M.vs.P.LumA.rs$rs$fdr < 0.01]
# cell.type.vec[M.vs.P.LumA.rs$rs$fdr < 0.05]
# 


save(file = 'client-side/output/TME.between.primary.and.metastatic.R.output/TME.between.primary.and.metastatic.RData',
     list=c('M.vs.P.Basal.rs','M.vs.P.Her2.rs','M.vs.P.LumB.rs','immune.score.matrix')
)


gene.symbol <- c('CD8A','CD8B','GZMA','GZMB','PRF1')
cell.type   <- 'CD8.T.cell'
ensemble.id <- c('ENSG00000153563','ENSG00000172116','ENSG00000145649','ENSG00000100453','ENSG00000180644')
CD8.peng.jiang.signature <- data.frame(gene.symbol = gene.symbol,cell.type=cell.type,ensemble.id = ensemble.id)


cell.type.vec        <- CD8.peng.jiang.signature$cell.type %>% unique
immune.score.matrix  <- compute.immune.score(CD8.peng.jiang.signature)

MET500.sample        <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample          <- pure.TCGA.breast.cancer.polyA.Basal.sample
M.vs.P.Basal.rs      <- compare.TME(.case.sample = MET500.sample,.control.sample = TCGA.sample,.immune.score.matrix = immune.score.matrix)
cell.type.vec[M.vs.P.Basal.rs$rs$fdr < 0.01]
cell.type.vec[M.vs.P.Basal.rs$rs$fdr < 0.05]


MET500.sample        <- MET500.breast.cancer.polyA.LumB.sample
MET500.sample        <- intersect(MET500.sample,MET500.liver.sample)
MET500.sample        <- MET500.sample[hepatocyte.abundance.vec[MET500.sample] < 0.2]
TCGA.sample          <- pure.TCGA.breast.cancer.polyA.LumB.sample
M.vs.P.LumB.rs       <- compare.TME(.case.sample = MET500.sample,.control.sample = TCGA.sample,.immune.score.matrix = immune.score.matrix)
cell.type.vec[M.vs.P.LumB.rs$rs$fdr < 0.01]
cell.type.vec[M.vs.P.LumB.rs$rs$fdr < 0.05]

M.vs.P.Basal.rs.with.peng.jiang.signature <- M.vs.P.Basal.rs
M.vs.P.LumB.rs.with.peng.jiang.signature  <- M.vs.P.LumB.rs


save(file = 'client-side/output/TME.between.primary.and.metastatic.R.output/TME.between.primary.and.metastatic.with.peng.jiang.T.cell.signature.RData',
     list=c('M.vs.P.Basal.rs.with.peng.jiang.signature','M.vs.P.LumB.rs.with.peng.jiang.signature','immune.score.matrix')
)
# ################# TME comaprision among different subtypes  ########################
# 
# 
# ####### well. for Basal.vs.other subtype, we need to use breast.basal.subtype.immune.cell.signature signuature, 
# ####### in this signature we removed genes which show DE between Basal and Luminal cell types
# immune.score.matrix  <- compute.immune.score(immune.cell.signature)
# MET500.liver.sample  <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# x                    <- hepatocyte.abundance.vec[MET500.liver.sample]
# x                    <- x[is.na(x) == FALSE]
# MET500.liver.sample <- names(x)[x < 0.2]
# 
# P.Basal.vs.LumB.purity.scaled   <- compare.TME(.case.sample = pure.TCGA.breast.cancer.polyA.Basal.sample,.control.sample = pure.TCGA.breast.cancer.polyA.LumB.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE)
# P.Basal.vs.LumB.purity.unscaled <- compare.TME(.case.sample = pure.TCGA.breast.cancer.polyA.Basal.sample,.control.sample = pure.TCGA.breast.cancer.polyA.LumB.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE)
# plot(x= -1 * (P.Basal.vs.LumB.purity.unscaled$rs$fdr %>% log10),y=-1 * (P.Basal.vs.LumB.purity.scaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[P.Basal.vs.LumB.purity.unscaled$rs$fdr < 0.05 & P.Basal.vs.LumB.purity.scaled$rs$fdr < 0.05]
# 
# M.Basal.vs.LumB.purity.scaled   <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Basal.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumB.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE,fit.method = 'lm')
# M.Basal.vs.LumB.purity.unscaled <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Basal.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumB.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE,fit.method = 'lm')
# plot(x= -1 * (M.Basal.vs.LumB.purity.scaled$rs$fdr %>% log10),y=-1 * (M.Basal.vs.LumB.purity.unscaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[M.Basal.vs.LumB.purity.scaled$rs$fdr < 0.01 & M.Basal.vs.LumB.purity.unscaled$rs$fdr < 0.01]
# 
# 
# 
# P.Basal.vs.LumA.purity.scaled   <- compare.TME(.case.sample = pure.TCGA.breast.cancer.polyA.Basal.sample,.control.sample = pure.TCGA.breast.cancer.polyA.LumA.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE)
# P.Basal.vs.LumA.purity.unscaled <- compare.TME(.case.sample = pure.TCGA.breast.cancer.polyA.Basal.sample,.control.sample = pure.TCGA.breast.cancer.polyA.LumA.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE)
# plot(x= -1 * (P.Basal.vs.LumA.purity.unscaled$rs$fdr %>% log10),y=-1 * (P.Basal.vs.LumA.purity.scaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[P.Basal.vs.LumA.purity.unscaled$rs$fdr < 0.01 & P.Basal.vs.LumA.purity.scaled$rs$fdr < 0.01]
# 
# MET500.liver.sample             <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# M.Basal.vs.LumA.purity.scaled   <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Basal.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumA.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE,fit.method = 'lm')
# M.Basal.vs.LumA.purity.unscaled <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Basal.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumA.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE,fit.method = 'lm')
# plot(x= -1 * (M.Basal.vs.LumA.purity.scaled$rs$fdr %>% log10),y=-1 * (M.Basal.vs.LumA.purity.unscaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[M.Basal.vs.LumB.purity.scaled$rs$fdr < 0.01 & M.Basal.vs.LumB.purity.unscaled$rs$fdr < 0.01]
# 
# 
# 
# P.Basal.vs.Her2.purity.scaled   <- compare.TME(.case.sample = pure.TCGA.breast.cancer.polyA.Basal.sample,.control.sample = pure.TCGA.breast.cancer.polyA.Her2.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE)
# P.Basal.vs.Her2.purity.unscaled <- compare.TME(.case.sample = pure.TCGA.breast.cancer.polyA.Basal.sample,.control.sample = pure.TCGA.breast.cancer.polyA.Her2.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE)
# plot(x= -1 * (P.Basal.vs.Her2.purity.unscaled$rs$fdr %>% log10),y=-1 * (P.Basal.vs.Her2.purity.scaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[P.Basal.vs.Her2.purity.unscaled$rs$fdr < 0.01 & P.Basal.vs.Her2.purity.scaled$rs$fdr < 0.01]
# 
# 
# MET500.liver.sample             <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# M.Basal.vs.Her2.purity.scaled   <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Basal.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.Her2.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE,fit.method = 'lm')
# M.Basal.vs.Her2.purity.unscaled <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Basal.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.Her2.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE,fit.method = 'lm')
# plot(x= -1 * (M.Basal.vs.Her2.purity.scaled$rs$fdr %>% log10),y=-1 * (M.Basal.vs.Her2.purity.unscaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[M.Basal.vs.Her2.purity.scaled$rs$fdr < 0.01 & M.Basal.vs.Her2.purity.unscaled$rs$fdr < 0.01]
# 
# 
# 
# 
# 
# 
# ####### well. for non-Basal subtype comparision, we assume they have the same cell type of origin ##########
# ####### just use the signature immune.cell.signature ##################
# immune.score.matrix  <- compute.immune.score(immune.cell.signature)
# 
# 
# P.Her2.vs.LumB.purity.scaled   <- compare.TME(.case.sample = TCGA.breast.cancer.polyA.Her2.sample,.control.sample = TCGA.breast.cancer.polyA.LumB.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE)
# P.Her2.vs.LumB.purity.unscaled <- compare.TME(.case.sample = TCGA.breast.cancer.polyA.Her2.sample,.control.sample = TCGA.breast.cancer.polyA.LumB.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE)
# plot(x= -1 * (P.Her2.vs.LumB.purity.unscaled$rs$fdr %>% log10),y=-1 * (P.Her2.vs.LumB.purity.scaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[P.Her2.vs.LumB.purity.unscaled$rs$fdr < 0.01 & P.Her2.vs.LumB.purity.scaled$rs$fdr < 0.01]
# 
# MET500.liver.sample            <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# M.Her2.vs.LumB.purity.scaled   <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Her2.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumB.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE, fit.method = 'lm')
# M.Her2.vs.LumB.purity.unscaled <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Her2.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumB.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE,fit.method = 'lm')
# plot(x= -1 * (M.Her2.vs.LumB.purity.scaled$rs$fdr %>% log10),y=-1 * (M.Her2.vs.LumB.purity.unscaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[M.Her2.vs.LumB.purity.scaled$rs$fdr < 0.01 & M.Her2.vs.LumB.purity.unscaled$rs$fdr < 0.01]
# 
# 
# 
# P.Her2.vs.LumA.purity.scaled   <- compare.TME(.case.sample = TCGA.breast.cancer.polyA.Her2.sample,.control.sample = TCGA.breast.cancer.polyA.LumA.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE)
# P.Her2.vs.LumA.purity.unscaled <- compare.TME(.case.sample = TCGA.breast.cancer.polyA.Her2.sample,.control.sample = TCGA.breast.cancer.polyA.LumA.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE)
# plot(x= -1 * (P.Her2.vs.LumA.purity.unscaled$rs$fdr %>% log10),y=-1 * (P.Her2.vs.LumA.purity.scaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[P.Her2.vs.LumA.purity.unscaled$rs$fdr < 0.01 & P.Her2.vs.LumA.purity.scaled$rs$fdr < 0.01]
# 
# MET500.liver.sample            <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# M.Her2.vs.LumA.purity.scaled   <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Her2.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumA.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE, fit.method = 'lm')
# M.Her2.vs.LumA.purity.unscaled <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.Her2.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumA.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE,fit.method = 'lm')
# plot(x= -1 * (M.Her2.vs.LumA.purity.scaled$rs$fdr %>% log10),y=-1 * (M.Her2.vs.LumA.purity.unscaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[M.Her2.vs.LumA.purity.scaled$rs$fdr < 0.01 & M.Her2.vs.LumA.purity.unscaled$rs$fdr < 0.01]
# 
# 
# 
# P.LumB.vs.LumA.purity.scaled   <- compare.TME(.case.sample = TCGA.breast.cancer.polyA.LumB.sample,.control.sample = TCGA.breast.cancer.polyA.LumA.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE)
# P.LumB.vs.LumA.purity.unscaled <- compare.TME(.case.sample = TCGA.breast.cancer.polyA.LumB.sample,.control.sample = TCGA.breast.cancer.polyA.LumA.sample,.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE)
# plot(x= -1 * (P.LumB.vs.LumA.purity.unscaled$rs$fdr %>% log10),y=-1 * (P.LumB.vs.LumA.purity.scaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[P.LumB.vs.LumA.purity.unscaled$rs$fdr < 0.01 & P.LumB.vs.LumA.purity.scaled$rs$fdr < 0.01]
# 
# MET500.liver.sample            <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# M.LumB.vs.LumA.purity.scaled   <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.LumB.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumA.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = TRUE, fit.method = 'lm')
# M.LumB.vs.LumA.purity.unscaled <- compare.TME(.case.sample = intersect(MET500.breast.cancer.polyA.LumB.sample,MET500.liver.sample),.control.sample = intersect(MET500.breast.cancer.polyA.LumA.sample,MET500.liver.sample),.immune.score.matrix = immune.score.matrix,.purity.scale = FALSE,fit.method = 'lm')
# plot(x= -1 * (M.LumB.vs.LumA.purity.scaled$rs$fdr %>% log10),y=-1 * (M.LumB.vs.LumA.purity.unscaled$rs$fdr %>% log10),xlab='fdr.unscaled',ylab='fdr.scaled')
# abline(v=2)
# abline(h=2)
# cell.type.vec[M.LumB.vs.LumA.purity.scaled$rs$fdr < 0.01 & M.LumB.vs.LumA.purity.unscaled$rs$fdr < 0.01]


# save(file = 'client-side/output/TME.between.primary.and.metastatic.R.output/TME.between.primary.and.metastatic.RData',
#      list=c('M.vs.P.Basal.rs','M.vs.P.Her2.rs','M.vs.P.LumA.rs','M.vs.P.LumB.rs',
#             'P.Basal.vs.Her2.purity.scaled',  'P.Basal.vs.LumA.purity.scaled',  'P.Basal.vs.LumB.purity.scaled',  'P.Her2.vs.LumA.purity.scaled',  'P.Her2.vs.LumB.purity.scaled',  'P.LumB.vs.LumA.purity.scaled',
#             'P.Basal.vs.Her2.purity.unscaled','P.Basal.vs.LumA.purity.unscaled','P.Basal.vs.LumB.purity.unscaled','P.Her2.vs.LumA.purity.unscaled','P.Her2.vs.LumB.purity.unscaled','P.LumB.vs.LumA.purity.unscaled',
#             'M.Basal.vs.Her2.purity.scaled',  'M.Basal.vs.LumA.purity.scaled',  'M.Basal.vs.LumB.purity.scaled',  'M.Her2.vs.LumA.purity.scaled',  'M.Her2.vs.LumB.purity.scaled',  'M.LumB.vs.LumA.purity.scaled',
#             'M.Basal.vs.Her2.purity.unscaled','M.Basal.vs.LumA.purity.unscaled','M.Basal.vs.LumB.purity.unscaled','M.Her2.vs.LumA.purity.unscaled','M.Her2.vs.LumB.purity.unscaled','M.LumB.vs.LumA.purity.unscaled'
#            )
#     )




################################################ Trash below this line ########################################################################################



################################################################################################
################ 
################ Compare TME difference between different subtypes, same site ##############
################ 
################################################################################################

# ############ check the marker gene expression in single-cell RNASeq dataset,they are NOT DE between different cell types ################
# load('server-side/RData/GSE113197_FPKM.RData')
# 
# apply(GSE113197[marker.gene,],2,median)
# 
# GSE113197.matrix <- as.matrix(GSE113197)
# 
# x <- GSE113197.matrix['ENSG00000186847',]+1  
# 
# 
# apply(GSE113197.matrix[marker.gene,],1,median)
# 
# GSE113197.matrix[marker.gene,] %>% t %>% boxplot
# 
# CD14 <- 'ENSG00000170458'
# KRT14 <- 'ENSG00000186847'
# plot(GSE113197.matrix[c(CD14,KRT14),] %>% t)
# 
# 
# ####### Maybe I could try CPE here, from Dvir's computation!!!
# ####### TME comparision between subtypes, primary site###########
# load("/Users/liuke/Project/InSilicoCRISPR/client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData")
# 
# TCGA.sample                      <- intersect(TCGA.breast.cancer.polyA.LumB.sample,TCGA.CD8.T.cell.level.df %>% rownames)
# # rs.TCGA             <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # purity.TCGA.LumB    <- rs.TCGA$correlation.matrix[,'BT483_BREAST']
# # purity.TCGA.LumB    <- scale(purity.TCGA.LumB,center = TRUE,scale = TRUE) %>% c
# purity.TCGA.LumB.CPE        <- TCGA.CD8.T.cell.level.df[TCGA.sample,'tumor.purity']
# names(purity.TCGA.LumB.CPE) <- TCGA.sample
# 
# TCGA.sample                       <- intersect(TCGA.breast.cancer.polyA.Basal.sample,TCGA.CD8.T.cell.level.df %>% rownames)
# # rs.TCGA             <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # purity.TCGA.Basal   <- rs.TCGA$correlation.matrix[,'HCC1569_BREAST']
# # purity.TCGA.Basal    <- scale(purity.TCGA.Basal,center = TRUE,scale = TRUE) %>% c
# purity.TCGA.Basal.CPE        <- TCGA.CD8.T.cell.level.df[TCGA.sample,'tumor.purity']
# names(purity.TCGA.Basal.CPE) <- TCGA.sample
# 
# 
# 
# 
# TCGA.sample         <- TCGA.breast.cancer.polyA.Her2.sample
# # rs.TCGA             <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # purity.TCGA.Her2    <- rs.TCGA$correlation.matrix[,'EFM192A_BREAST']
# # purity.TCGA.Her2.CPE        <- TCGA.CD8.T.cell.level.df[TCGA.sample,'tumor.purity']
# purity.TCGA.Her2.CPE        <- TCGA.CD8.T.cell.level.df[TCGA.sample,'tumor.purity']
# names(purity.TCGA.Her2.CPE) <- TCGA.sample
# 
# 
# TCGA.sample         <- TCGA.breast.cancer.polyA.LumA.sample
# # rs.TCGA             <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # purity.TCGA.Her2    <- rs.TCGA$correlation.matrix[,'EFM192A_BREAST']
# # purity.TCGA.Her2.CPE        <- TCGA.CD8.T.cell.level.df[TCGA.sample,'tumor.purity']
# purity.TCGA.LumA.CPE        <- TCGA.CD8.T.cell.level.df[TCGA.sample,'tumor.purity']
# names(purity.TCGA.LumA.CPE) <- TCGA.sample
# 
# 
# 
# 
# x             <- c(purity.TCGA.LumB.CPE,purity.TCGA.Basal.CPE)
# rs <- foreach(g= marker.gene,.combine='rbind') %do% {
#     df                <- data.frame(x=x,y=cancer.data[g,c(names(purity.TCGA.LumB.CPE),names(purity.TCGA.Basal.CPE))])
#     df                <- df[complete.cases(df),]
#     df$subtype        <- ifelse(rownames(df) %in% TCGA.breast.cancer.polyA.LumB.sample,'LumB','Basal')
#     df$subtype        <- relevel(df$subtype %>% factor,ref = 'LumB')
#     fit.rs            <- lm(df,formula = y ~ x + subtype)
#     r                 <- summary(fit.rs)$coefficients[3,c(1,4)]
#     
#     # loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#     # adjusted.expr     <- loess.fit$residuals
#     # p.value           <- wilcox.test(adjusted.expr[names(purity.TCGA.Basal.CPE)],adjusted.expr[names(purity.TCGA.LumB.CPE)])$p.value
#     # delta             <- median(adjusted.expr[names(purity.TCGA.Basal.CPE)]) -  median(adjusted.expr[names(purity.TCGA.LumB.CPE)])
#     data.frame(effect.size=r[1],cancer.p.value=r[2])
# }
# rownames(rs)          <- marker.gene.symbol
# Basal.vs.LumB.rs      <- rs
# 
# 
# g <- 'ENSG00000170458'   #CD14
# g <- 'ENSG00000049768'   #CD68
# x             <- c(purity.TCGA.LumB.CPE,purity.TCGA.Basal.CPE)
# df            <- data.frame(x=x,y=cancer.data[g,names(x)])
# df$data.source <- ifelse(rownames(df) %in% TCGA.breast.cancer.polyA.LumB.sample,'LumB','Basal')
# ggplot(df,aes(x=x,y=y,color=factor(data.source))) + geom_point(aes(color=factor(data.source)),size=4) + geom_smooth(method = 'lm') + ggplot.style + xlab('purity') + ylab('gene expression') 
# 
# 
# 
# 
# 
# x             <- c(purity.TCGA.Her2.CPE,purity.TCGA.Basal.CPE)
# rs <- foreach(g= marker.gene,.combine='rbind') %do% {
#   df                <- data.frame(x=x,y=cancer.data[g,c(names(purity.TCGA.Her2.CPE),names(purity.TCGA.Basal.CPE))])
#   df                <- df[complete.cases(df),]
#   df$subtype        <- ifelse(rownames(df) %in% TCGA.breast.cancer.polyA.Her2.sample,'Her2','Basal')
#   df$subtype        <- relevel(df$subtype %>% factor,ref = 'Her2')
#   fit.rs            <- lm(df,formula = y ~ x + subtype)
#   r                 <- summary(fit.rs)$coefficients[3,c(1,4)]
#   
#   # loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#   # adjusted.expr     <- loess.fit$residuals
#   # p.value           <- wilcox.test(adjusted.expr[names(purity.TCGA.Basal.CPE)],adjusted.expr[names(purity.TCGA.LumB.CPE)])$p.value
#   # delta             <- median(adjusted.expr[names(purity.TCGA.Basal.CPE)]) -  median(adjusted.expr[names(purity.TCGA.LumB.CPE)])
#   data.frame(effect.size=r[1],cancer.p.value=r[2])
# }
# rownames(rs)          <- marker.gene.symbol
# Basal.vs.Her2.rs      <- rs
# 
# 
# g <- 'ENSG00000170458'   #CD14
# x             <- c(purity.TCGA.Her2.CPE,purity.TCGA.Basal.CPE)
# df            <- data.frame(x=x,y=cancer.data[g,names(x)])
# df$data.source <- ifelse(rownames(df) %in% TCGA.breast.cancer.polyA.Her2.sample,'Her2','Basal')
# ggplot(df,aes(x=x,y=y)) + geom_point(aes(color=factor(data.source)),size=4) + geom_smooth(method = 'lm') + ggplot.style + xlab('purity') + ylab('gene expression') +  theme(legend.position = "none")
# 
# 
# 
# ################################# MET500, TME comparision between different subtypes ##############
# 
# 
# MET500.sample <- MET500.breast.cancer.polyA.Basal.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500.Basal <- rs.MET500$correlation.matrix[,'HCC70_BREAST'] # use HCC70 cell line
# 
# 
# MET500.sample <- MET500.breast.cancer.polyA.LumB.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500.LumB <- rs.MET500$correlation.matrix[,'BT483_BREAST'] # use BT483_BREAST cell line
# 
# 
# 
# MET500.sample <- MET500.breast.cancer.polyA.LumA.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500.LumA <- rs.MET500$correlation.matrix[,'EFM192A_BREAST'] # use EFM192A_BREAST cell line
# 
# 
# 
# 
# 
# x             <- c(purity.MET500.Basal,purity.MET500.LumB,purity.MET500.LumA)
# rs <- foreach(g= marker.gene,.combine='rbind') %do% {
#   df                <- data.frame(x=x,y=MET500.log2.fpkm.matrix[g,names(x)])
#   df                <- df[complete.cases(df),]
#   df$subtype        <- ifelse(rownames(df) %in% MET500.breast.cancer.polyA.Basal.sample,'Basal','Other')
#   df$subtype        <- relevel(df$subtype %>% factor,ref = 'Other')
#   fit.rs            <- lm(df,formula = y ~ x + subtype)
#   r                 <- summary(fit.rs)$coefficients[3,c(1,4)]
#   
#   # loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#   # adjusted.expr     <- loess.fit$residuals
#   # p.value           <- wilcox.test(adjusted.expr[names(purity.TCGA.Basal.CPE)],adjusted.expr[names(purity.TCGA.LumB.CPE)])$p.value
#   # delta             <- median(adjusted.expr[names(purity.TCGA.Basal.CPE)]) -  median(adjusted.expr[names(purity.TCGA.LumB.CPE)])
#   data.frame(effect.size=r[1],cancer.p.value=r[2])
# }
# rownames(rs)          <- marker.gene.symbol
# Basal.vs.LumB.rs      <- rs
# 
# 
# 
# g             <- 'ENSG00000170458'   #CD14
# x             <- c(purity.MET500.Basal,purity.MET500.LumB,purity.MET500.LumA)
# df            <- data.frame(x=x,y=MET500.log2.fpkm.matrix[g,names(x)])
# df$data.source <- ifelse(rownames(df) %in% MET500.breast.cancer.polyA.Basal.sample,'Basal','Other')
# ggplot(df,aes(x=x,y=y)) + geom_point(aes(color=factor(data.source)),size=4) + geom_smooth(method = 'lm') + ggplot.style + xlab('purity') + ylab('gene expression') +  theme(legend.position = "none")
# 
# 
# 
# 
# library(xCell)
# library(plyr)
# library(dplyr)
# library (VennDiagram)
# library(gplots)
# library(foreach)
# library(org.Hs.eg.db)
# library(pheatmap)
# library(foreach)
# library(AnnotationDbi)
# source('client-side/code/util.R')
# 
# #Ref: RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types
# 
# ######################### Comparsion between GTex Liver and Breast, remove the DE genes that may confound our further analysis ###################
# load('server-side//RData//Breast - Mammary Tissue.RData')
# idx                                <- sample.meta.df$sample.id[sample.meta.df$gender == 'Female'] %>% as.character
# flag                               <- colnames(log2.read.count.matrix) %in% idx
# GTex.breast.log2.fpkm.matrix       <- log2.fpkm.matrix
# 
# load('server-side//RData//Liver.RData')
# idx                                <- sample.meta.df$sample.id[sample.meta.df$gender == 'Female'] %>% as.character
# flag                               <- colnames(log2.read.count.matrix) %in% idx
# GTex.liver.log2.fpkm.matrix        <- log2.fpkm.matrix[,flag]
# 
# ######## map ENSEMBLE id to gene symbol ##############
# protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
# ensemble.gene.id       <- intersect(rownames(GTex.liver.log2.fpkm.matrix),protein.coding.gene.id)
# 
# df    <- AnnotationDbi::select(org.Hs.eg.db,keys=ensemble.gene.id, columns=c('SYMBOL'), keytype="ENSEMBL")
# df    <- df[complete.cases(df),]
# tmp   <- ddply(df,.(ENSEMBL),nrow)
# tmp   <- tmp[order(tmp$V1,decreasing = TRUE),]
# g1    <- tmp$ENSEMBL[tmp$V1 == 1] %>% as.character()
# tmp   <- ddply(df,.(SYMBOL),nrow)
# tmp   <- tmp[order(tmp$V1,decreasing = TRUE),]
# g2    <- tmp$SYMBOL[tmp$V1 == 1] %>% as.character()
# mapping.df           <- df[(df$ENSEMBL %in% g1)  & (df$SYMBOL %in% g2),]
# rownames(mapping.df) <- mapping.df$ENSEMBL
# 
# 
# 
# data                  <- cbind(GTex.breast.log2.fpkm.matrix,GTex.liver.log2.fpkm.matrix[rownames(GTex.breast.log2.fpkm.matrix),])
# data                  <- data[rownames(data) %in% mapping.df$ENSEMBL,]
# rownames(data)        <- mapping.df[rownames(data),'SYMBOL']
# GTex.data.for.xCell   <- data
# 
# GTex.data.for.xCell.rank      <- apply(GTex.data.for.xCell,2,rank)
# breast.rank.data              <- GTex.data.for.xCell.rank[,colnames(GTex.breast.log2.fpkm.matrix)]
# liver.rank.data               <- GTex.data.for.xCell.rank[,colnames(GTex.liver.log2.fpkm.matrix)]
# 
# rank.de.df <- foreach(i = 1:nrow(breast.rank.data),.combine = 'rbind') %do% {
#     data.frame(delta = median(liver.rank.data[i,]) - median(breast.rank.data[i,]),p.value=wilcox.test(breast.rank.data[i,],liver.rank.data[i,])$p.value) 
# }
# rownames(rank.de.df) <- rownames(breast.rank.data)
# rank.de.df$fdr       <- p.adjust(rank.de.df$p.value,method='fdr')
# flag                 <- (rank.de.df$fdr < 0.01) 
# de.gene              <- rownames(rank.de.df)[flag]
# 
# 
# xCell.corrected.gene.set <- foreach(x=xCell.data$signatures) %do% {
#     setdiff(x@geneIds,de.gene) 
# }
# names(xCell.corrected.gene.set) <- names(xCell.data$signatures) ############# remove DE genes from xCell signatures. Maybe here we should also remove genes DE between metastatic and primary 
# 
# L <- list()
# foreach(i = names(xCell.corrected.gene.set)) %do% {
#     tmp <- strsplit(x = i,split='%')[[1]] %>% unlist      
#     cell.type <- tmp[1]
#     L[[cell.type]] <- c(L[[cell.type]],xCell.corrected.gene.set[[i]])
#     L[[cell.type]] <- unique(L[[cell.type]])
# }
# xCell.corrected.gene.set <- L
# 
# 
# ## maybe here I should filter out the genes not expressed
# 
# 
# 
# GTex.ssgsea.scores <- GSVA::gsva(expr=GTex.data.for.xCell,
#                                  xCell.corrected.gene.set, method = "ssgsea",
#                                  ssgsea.norm = FALSE)
# breast.data <- GTex.ssgsea.scores[,colnames(GTex.breast.log2.fpkm.matrix)]
# liver.data  <- GTex.ssgsea.scores[,colnames(GTex.liver.log2.fpkm.matrix)]
# GTex.de.df  <- foreach(i = 1:nrow(GTex.ssgsea.scores),.combine='rbind') %do% {
#     data.frame(delta = median(liver.data[i,]) - median(breast.data[i,]),p.value= wilcox.test(breast.data[i,],liver.data[i,])$p.value)
# }
# rownames(GTex.de.df) <- rownames(breast.data)
# GTex.de.df$fdr       <- p.adjust(GTex.de.df$p.value,method='fdr') ## Run ssGSEA with corrected dataset, no significant geneset
# sum(GTex.de.df$fdr < 0.05)
# GTex.corrected.gene.set.da.df <- GTex.de.df
# 
# 
# GTex.ssgsea.scores <- GSVA::gsva(expr=GTex.data.for.xCell,
#                                  xCell.data$signatures, method = "ssgsea",
#                                  ssgsea.norm = FALSE)
# breast.data <- GTex.ssgsea.scores[,colnames(GTex.breast.log2.fpkm.matrix)]
# liver.data  <- GTex.ssgsea.scores[,colnames(GTex.liver.log2.fpkm.matrix)]
# GTex.de.df  <- foreach(i = 1:nrow(GTex.ssgsea.scores),.combine='rbind') %do% {
#     data.frame(delta = median(liver.data[i,]) - median(breast.data[i,]),p.value= wilcox.test(breast.data[i,],liver.data[i,])$p.value)
# }
# rownames(GTex.de.df) <- rownames(breast.data)
# GTex.de.df$fdr       <- p.adjust(GTex.de.df$p.value,method='fdr') ## Run ssGSEA with corrected dataset, no significant geneset
# sum(GTex.de.df$fdr < 0.05)
# GTex.gene.set.da.df  <- GTex.de.df
# 
# 
# boxplot(GTex.gene.set.da.df$p.value,GTex.corrected.gene.set.da.df$p.value)
# hist(GTex.gene.set.da.df$p.value)
# hist(GTex.corrected.gene.set.da.df$p.value)
# 
# 
# ############# compute correlation values with cell line, as estimation of tumor purity
# load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')
# CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
# CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
# tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
# tmp.rank                    <- apply(tmp,2,rank)
# rank.mean                   <- apply(tmp.rank,1,mean)
# rank.sd                     <- apply(tmp.rank,1,sd)
# plot(x=rank.mean,y=rank.sd)
# lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
# CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
# 
# 
# 
# 
# ####################### prepare cancer data for xCell analysis ###########################
# load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
# load('server-side/RData//Breast Invasive Carcinoma.RData')
# load('~/Project/Cancer2CellLine/client-side/output//MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
# load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
# 
# TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
# 
# data                  <- cbind(TCGA.breast.cancer.log2.fpkm.matrix,MET500.log2.fpkm.matrix[rownames(TCGA.breast.cancer.log2.fpkm.matrix),])
# data                  <- data[rownames(data) %in% mapping.df$ENSEMBL,]
# rownames(data)        <- mapping.df[rownames(data),'SYMBOL']
# cancer.data.for.xCell <- data
# 
# MET500.liver.sample   <- rownames(MET500.sample.meta)[MET500.sample.meta$biopsy.site == 'LIVER'] 
# 
# 
# ################# Basal subtype  ########################
# MET500.sample <- MET500.breast.cancer.polyA.Basal.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# TCGA.sample   <- TCGA.breast.cancer.polyA.Basal.sample
# 
# cancer.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],
#                                    xCell.corrected.gene.set, method = "ssgsea",
#                                    ssgsea.norm = FALSE) # hmm, maybe here should set some cutoff, ssgsea score smaller than this cut-off truncate to zero. 
# 
# #xCell.score.matrix <- xCellAnalysis(cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],file.name = 'client-side/output/TME.R.output/basal.raw.score.txt' )
# 
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500 <- rs.MET500$correlation.matrix[,'HCC70_BREAST'] # use HCC70 cell line
# purity.TCGA   <- rs.TCGA$correlation.matrix[,'HCC70_BREAST']
# x             <- c(purity.TCGA,purity.MET500)
# 
# 
# cancer.geneset.da.df <- foreach(i =1:nrow(cancer.ssgsea.scores),.combine='rbind') %do% {
#     df                <- data.frame(x=x,y=cancer.ssgsea.scores[i,])
#     loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#     adjusted.expr     <- loess.fit$residuals
#     p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
#     delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
#     data.frame(cancer.delta=delta,cancer.p.value=p.value)
#     
# }
# rownames(cancer.geneset.da.df)  <- rownames(cancer.ssgsea.scores)
# cancer.geneset.da.df$cancer.fdr <- p.adjust(cancer.geneset.da.df$cancer.p.value,method='fdr')
# da.gene.set                     <- rownames(cancer.geneset.da.df)[cancer.geneset.da.df$cancer.fdr < 0.05]
# 
# Basal.cancer.geneset.da.df      <- cancer.geneset.da.df
# Basal.da.gene.set               <- da.gene.set
# 
# plot(x=Basal.cancer.geneset.da.df$cancer.delta,y = -1 * log10(Basal.cancer.geneset.da.df$cancer.fdr))
# 
# ################# LumB subtype ########################
# MET500.sample <- MET500.breast.cancer.polyA.LumB.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# TCGA.sample   <- TCGA.breast.cancer.polyA.LumB.sample
# 
# cancer.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],
#                                    xCell.corrected.gene.set, method = "ssgsea",
#                                    ssgsea.norm = FALSE)
# 
# 
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500 <- rs.MET500$correlation.matrix[,'BT483_BREAST'] # use BT483 cell line
# purity.TCGA   <- rs.TCGA$correlation.matrix[,'BT483_BREAST']
# x             <- c(purity.TCGA,purity.MET500)
# 
# 
# cancer.geneset.da.df <- foreach(i =1:nrow(cancer.ssgsea.scores),.combine='rbind') %do% {
#     df                <- data.frame(x=x,y=cancer.ssgsea.scores[i,])
#     loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#     adjusted.expr     <- loess.fit$residuals
#     p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
#     delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
#     data.frame(cancer.delta=delta,cancer.p.value=p.value)
#   
# }
# rownames(cancer.geneset.da.df)  <- rownames(cancer.ssgsea.scores)
# cancer.geneset.da.df$cancer.fdr <- p.adjust(cancer.geneset.da.df$cancer.p.value,method='fdr')
# da.gene.set                     <- rownames(cancer.geneset.da.df)[cancer.geneset.da.df$cancer.fdr < 0.05]
# 
# LumB.cancer.geneset.da.df      <- cancer.geneset.da.df
# LumB.da.gene.set               <- da.gene.set
# 
# 
# 
# 
# ################# Her2 subtype ########################
# MET500.sample <- MET500.breast.cancer.polyA.Her2.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# TCGA.sample   <- TCGA.breast.cancer.polyA.Her2.sample
# 
# cancer.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],
#                                    xCell.corrected.gene.set, method = "ssgsea",
#                                    ssgsea.norm = FALSE)
# 
# 
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500 <- rs.MET500$correlation.matrix[,'EFM192A_BREAST'] # use EFM192A cell line
# purity.TCGA   <- rs.TCGA$correlation.matrix[,'EFM192A_BREAST']
# x             <- c(purity.TCGA,purity.MET500)
# 
# 
# cancer.geneset.da.df <- foreach(i =1:nrow(cancer.ssgsea.scores),.combine='rbind') %do% {
#     df                <- data.frame(x=x,y=cancer.ssgsea.scores[i,])
#     loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#     adjusted.expr     <- loess.fit$residuals
#     p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
#     delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
#     data.frame(cancer.delta=delta,cancer.p.value=p.value)
#   
# }
# rownames(cancer.geneset.da.df)  <- rownames(cancer.ssgsea.scores)
# cancer.geneset.da.df$cancer.fdr <- p.adjust(cancer.geneset.da.df$cancer.p.value,method='fdr')
# da.gene.set                     <- rownames(cancer.geneset.da.df)[cancer.geneset.da.df$cancer.fdr < 0.05]
# 
# Her2.cancer.geneset.da.df      <- cancer.geneset.da.df
# Her2.da.gene.set               <- da.gene.set
# 
# 
# 
# ####### Here, generate some random genesets to evaluate how large the batch effect is ###############
# m <- foreach( i = 1:10,.combine='cbind') %do% {
#     shuffle.gene.set <- foreach(x=xCell.corrected.gene.set) %do% {
#         sample(rownames(cancer.data.for.xCell),length(x)) 
#     }
#     names(shuffle.gene.set) <- names(xCell.corrected.gene.set)
# 
#     shuffle.cancer.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],
#                                            shuffle.gene.set[h], method = "ssgsea",
#                                            ssgsea.norm = FALSE)
#     x             <- c(purity.TCGA,purity.MET500)
# 
# shuffle.cancer.de.df <- foreach(i =1:nrow(shuffle.cancer.ssgsea.scores),.combine='rbind') %do% {
#     df                <- data.frame(x=x,y=shuffle.cancer.ssgsea.scores[i,])
#     loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#     adjusted.expr     <- loess.fit$residuals
#     p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
#     delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
#     data.frame(cancer.delta=delta,cancer.p.value=p.value)
# }
# #rownames(shuffle.cancer.de.df)  <- rownames(cancer.ssgsea.scores)
# shuffle.cancer.de.df$cancer.delta
# }
# 
# 
# 
# ############################ Now, let us try marker genes to analyze TME ##############################
# 
# marker.gene <- c('CD3D','CD4','CD8A','CD8B',    # T cell
#                  'CD19','MS4A1',                # CD19,CD20, B cell
#                  'NCAM1',                       # CD56, NK cell
#                  'ITGA2B'                       # CD41, Platlet
#                  )
# 
# 
# ##################################### Basal subtype #####################################
# MET500.sample <- MET500.breast.cancer.polyA.Basal.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# TCGA.sample   <- TCGA.breast.cancer.polyA.Basal.sample
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500 <- rs.MET500$correlation.matrix[,'HCC70_BREAST'] # use HCC70 cell line
# purity.TCGA   <- rs.TCGA$correlation.matrix[,'HCC70_BREAST']
# x             <- c(purity.TCGA,purity.MET500)
# 
# rs <- foreach(g= marker.gene,.combine='rbind') %do% {
#     df                <- data.frame(x=x,y=cancer.data.for.xCell[g,c(TCGA.sample,MET500.sample)])
#     loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#     adjusted.expr     <- loess.fit$residuals
#     p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
#     delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
#     data.frame(effect.size=delta,cancer.p.value=p.value)
# }
# rownames(rs) <- marker.gene
# Basal.rs <- rs
# 
# 
# ##################################### LumB subtype #####################################
# MET500.sample <- MET500.breast.cancer.polyA.LumB.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# TCGA.sample   <- TCGA.breast.cancer.polyA.LumB.sample
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500 <- rs.MET500$correlation.matrix[,'BT483_BREAST'] # use HCC70 cell line
# purity.TCGA   <- rs.TCGA$correlation.matrix[,'BT483_BREAST']
# x             <- c(purity.TCGA,purity.MET500)
# 
# rs <- foreach(g= marker.gene,.combine='rbind') %do% {
#     df                <- data.frame(x=x,y=cancer.data.for.xCell[g,c(TCGA.sample,MET500.sample)])
#     loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#     adjusted.expr     <- loess.fit$residuals
#     p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
#     delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
#     data.frame(effect.size=delta,cancer.p.value=p.value)
# }
# rownames(rs) <- marker.gene
# LumB.rs <- rs
# 
# 
# ##################################### Her2 subtype #####################################
# MET500.sample <- MET500.breast.cancer.polyA.Her2.sample
# MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# TCGA.sample   <- TCGA.breast.cancer.polyA.Her2.sample
# rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# purity.MET500 <- rs.MET500$correlation.matrix[,'EFM192A_BREAST'] # use HCC70 cell line
# purity.TCGA   <- rs.TCGA$correlation.matrix[,'EFM192A_BREAST']
# x             <- c(purity.TCGA,purity.MET500)
# 
# rs <- foreach(g= marker.gene,.combine='rbind') %do% {
#   df                <- data.frame(x=x,y=cancer.data.for.xCell[g,c(TCGA.sample,MET500.sample)])
#   loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
#   adjusted.expr     <- loess.fit$residuals
#   p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
#   delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
#   data.frame(effect.size=delta,cancer.p.value=p.value)
# }
# rownames(rs) <- marker.gene
# Her2.rs <- rs
# ##################################### Trash ######################################
# # primary.df   <- foreach(i = 1:nrow(cancer.ssgsea.scores),.combine='rbind') %do% {
# #   data.frame(p.value = wilcox.test(cancer.ssgsea.scores[i,TCGA.sample],breast.data[i,])$p.value,
# #              delta   = median(breast.data[i,]) -  median(cancer.ssgsea.scores[i,TCGA.sample])
# #             )
# # }
# # rownames(primary.df) <- rownames(cancer.ssgsea.scores)
# # 
# # metastatic.df <- foreach(i = 1:nrow(cancer.ssgsea.scores),.combine='rbind') %do% {
# #   data.frame(p.value=wilcox.test(cancer.ssgsea.scores[i,MET500.sample],liver.data[i,])$p.value,
# #              delta = median(liver.data[i,]) -  median(cancer.ssgsea.scores[i,MET500.sample])
# #   )
# #   
# # }
# # rownames(metastatic.df) <- rownames(cancer.ssgsea.scores)
# 
# # tmp <- cbind(primary.df,metastatic.df)
# # colnames(tmp) <- c('primary.p.value','primary.delta','metastatic.p.value','metastatic.delta')
# # 
# # tmp$primary.fdr    <- p.adjust(tmp$primary.p.value,method='fdr')
# # tmp$metastatic.fdr <- p.adjust(tmp$metastatic.p.value,method='fdr')
# # tmp$shape                                                       <- 'NS'
# # tmp[tmp$primary.fdr    < 0.05,'shape']                          <- 'PS'
# # tmp[tmp$metastatic.fdr < 0.05,'shape']                          <- 'MS'
# # tmp[tmp$primary.fdr < 0.05 & tmp$metastatic.fdr < 0.05,'shape'] <- 'BS'
# # 
# # 
# # ggplot(tmp,aes(x=primary.delta,y=metastatic.delta,color=factor(shape))) + geom_point()
# # 
# # 
# # 
# # 
# # plot(x=primary.p.value ,y=metastatic.p.value,xlim=c(-8000,8000),ylim=c(-8000,8000) )
# # lines(c(-8000,8000),c(-8000,8000))
# 
# # plot(x=GTex.liver.ssgsea.score.median - GTex.breast.ssgsea.score.median,y=p.value.vec)
# # plot(y=ssgsea.scores['CD8+ T-cells%BLUEPRINT%1.txt',],x=c(purity.TCGA,purity.MET500),col=ifelse(names(x) %in% MET500.sample,'red','black'),pch=19)
# # plot(y=ssgsea.scores['Astrocytes%ENCODE%1.txt',],x=c(purity.TCGA,purity.MET500),col=ifelse(names(x) %in% MET500.sample,'red','black'),pch=19)
# # 
# # plot(y=cancer.ssgsea.scores[i,],x=c(purity.TCGA,purity.MET500),col=ifelse(names(x) %in% MET500.sample,'red','black'),pch=19)
# # 
# # 
# # 
# # plot(y=cancer.ssgsea.scores['Class-switched memory B-cells',],x=c(purity.TCGA,purity.MET500),col=ifelse(names(x) %in% MET500.sample,'red','black'),pch=19)
# # i <- 'CD8+ T-cells'
# # plot(y=cancer.ssgsea.scores[i,],x=c(purity.TCGA,purity.MET500),col=ifelse(names(x) %in% MET500.sample,'red','black'),pch=19)
# # 
# # 
# # GZMB <- 'ENSG00000100453'
# # 
# # df                <- data.frame(x=x,y=cancer.data.for.xCell['CD3',c(TCGA.sample,MET500.sample)])
# # loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# # adjusted.expr     <- loess.fit$residuals
# # p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
# # delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
# # data.frame(cancer.delta=delta,cancer.p.value=p.value)
# # plot(x=df$x,y=df$y,col=ifelse(rownames(df) %in% MET500.sample,'red','black'),pch=19)
# # 
# # df                <- data.frame(x=x,y=cancer.data.for.xCell['CD4',c(TCGA.sample,MET500.sample)])
# # loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# # adjusted.expr     <- loess.fit$residuals
# # p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
# # delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
# # data.frame(cancer.delta=delta,cancer.p.value=p.value)
# # plot(x=df$x,y=df$y,col=ifelse(rownames(df) %in% MET500.sample,'red','black'),pch=19)
# # 
# # df                <- data.frame(x=x,y=cancer.data.for.xCell['CD8A',c(TCGA.sample,MET500.sample)])
# # loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# # adjusted.expr     <- loess.fit$residuals
# # p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
# # delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
# # data.frame(cancer.delta=delta,cancer.p.value=p.value)
# # plot(x=df$x,y=df$y,col=ifelse(rownames(df) %in% MET500.sample,'red','black'),pch=19)
# # 
# # 
# # ## NK Cell, CD56
# # df                <- data.frame(x=x,y=cancer.data.for.xCell['NCAM1',c(TCGA.sample,MET500.sample)])
# # loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# # adjusted.expr     <- loess.fit$residuals
# # p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
# # delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
# # data.frame(cancer.delta=delta,cancer.p.value=p.value)
# # plot(x=df$x,y=df$y,col=ifelse(rownames(df) %in% MET500.sample,'red','black'),pch=19)
# # 
# # ## Macrophage 
# # 
# # 
# # ## Platelet  CD41 (aka ITGA2B)
# # df                <- data.frame(x=x,y=cancer.data.for.xCell['ITGA2B',c(TCGA.sample,MET500.sample)])
# # loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# # adjusted.expr     <- loess.fit$residuals
# # p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
# # delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
# # data.frame(cancer.delta=delta,cancer.p.value=p.value)
# # plot(x=df$x,y=df$y,col=ifelse(rownames(df) %in% MET500.sample,'red','black'),pch=19)
# 
# 
# # ################# LumB subtype  ########################
# # MET500.sample <- MET500.breast.cancer.polyA.LumB.sample
# # MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# # TCGA.sample   <- TCGA.breast.cancer.polyA.LumB.sample
# # 
# # rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # purity.MET500 <- rs.MET500$correlation.matrix[,'BT483_BREAST'] # use BT483_BREAST cell line
# # purity.TCGA   <- rs.TCGA$correlation.matrix[,'BT483_BREAST']
# # x             <- c(purity.TCGA,purity.MET500)
# # 
# # 
# # rs <- foreach(g= marker.gene,.combine='rbind') %do% {
# #   df                <- data.frame(x=x,y=cancer.data[g,c(TCGA.sample,MET500.sample)])
# #   loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# #   adjusted.expr     <- loess.fit$residuals
# #   p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
# #   delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
# #   data.frame(effect.size=delta,cancer.p.value=p.value)
# # }
# # rownames(rs) <- marker.gene.symbol
# # LumB.rs      <- rs
# # LumB.rs$fdr  <- p.adjust(LumB.rs$cancer.p.value,method='fdr')
# # 
# # LumB.rs[LumB.rs$cancer.p.value < 0.05,]
# # 
# # 
# # g                 <- 'ENSG00000149294' # CD56, NK cell marker
# # 
# # df                <- data.frame(x=x,y=cancer.data[g,c(TCGA.sample,MET500.sample)])
# # df$data.source <- ifelse(rownames(df) %in% TCGA.sample,'TCGA','MET500')
# # ggplot(df,aes(x=x,y=y)) + geom_point(aes(color=factor(data.source)),size=4) + geom_smooth(method='loess',se=FALSE) + ggplot.style +  theme(legend.position = "none") + xlab('purity') + ylab('gene expression')
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ################# Her2 subtype  ########################
# # MET500.sample <- MET500.breast.cancer.polyA.Her2.sample
# # MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# # TCGA.sample   <- TCGA.breast.cancer.polyA.Her2.sample
# # 
# # rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # purity.MET500 <- rs.MET500$correlation.matrix[,'EFM192A_BREAST'] # use EFM192A_BREAST cell line
# # purity.TCGA   <- rs.TCGA$correlation.matrix[,'EFM192A_BREAST']
# # x             <- c(purity.TCGA,purity.MET500)
# # 
# # 
# # rs <- foreach(g= marker.gene,.combine='rbind') %do% {
# #   df                <- data.frame(x=x,y=cancer.data[g,c(TCGA.sample,MET500.sample)])
# #   loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# #   adjusted.expr     <- loess.fit$residuals
# #   p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
# #   delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
# #   data.frame(effect.size=delta,cancer.p.value=p.value)
# # }
# # rownames(rs) <- marker.gene.symbol
# # Her2.rs      <- rs
# # Her2.rs$fdr  <- p.adjust(Her2.rs$cancer.p.value,method='fdr')
# # 
# # 
# # 
# # 
# # Her2.rs[Her2.rs$cancer.p.value < 0.05,]
# # 
# # g <- 'ENSG00000167286' # CD3D
# # df                <- data.frame(x=x,y=cancer.data[g,c(TCGA.sample,MET500.sample)])
# # df$data.source <- ifelse(rownames(df) %in% TCGA.sample,'TCGA','MET500')
# # ggplot(df,aes(x=x,y=y)) + geom_point(aes(color=factor(data.source)),size=4) + geom_smooth(method='loess',se=FALSE) + ggplot.style +  theme(legend.position = "none") + xlab('purity') + ylab('gene expression')
# # 
# # ################# LumA subtype  ########################
# # MET500.sample <- MET500.breast.cancer.polyA.LumA.sample
# # MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
# # TCGA.sample   <- TCGA.breast.cancer.polyA.LumA.sample
# # 
# # rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
# # purity.MET500 <- rs.MET500$correlation.matrix[,'MDAMB415_BREAST'] # use MDAMB415_BREAST cell line
# # purity.TCGA   <- rs.TCGA$correlation.matrix[,'MDAMB415_BREAST']
# # x             <- c(purity.TCGA,purity.MET500)
# # 
# # 
# # rs <- foreach(g= marker.gene,.combine='rbind') %do% {
# #   df                <- data.frame(x=x,y=cancer.data[g,c(TCGA.sample,MET500.sample)])
# #   loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
# #   adjusted.expr     <- loess.fit$residuals
# #   p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
# #   delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
# #   data.frame(effect.size=delta,cancer.p.value=p.value)
# # }
# # rownames(rs) <- marker.gene.symbol
# # LumA.rs      <- rs
# # LumA.rs$fdr  <- p.adjust(LumA.rs$cancer.p.value,method='fdr')
# # 
# # LumA.rs[LumA.rs$cancer.p.value < 0.05,]
# # 
# # 
# # g                 <- 'ENSG00000172116' #CD8B
# # df                <- data.frame(x=x,y=cancer.data[g,c(TCGA.sample,MET500.sample)])
# # df$data.source <- ifelse(rownames(df) %in% TCGA.sample,'TCGA','MET500')
# # ggplot(df,aes(x=x,y=y)) + geom_point(aes(color=factor(data.source)),size=4) + geom_smooth(method='loess',se=FALSE) + ggplot.style +  theme(legend.position = "none") + xlab('purity') + ylab('gene expression')
# # 
# # 
# # 
# # combined.p.value.vec <- foreach(g= marker.gene.symbol,.combine = 'c') %do% {
# #   s <- -2 * log(LumA.rs[g,'cancer.p.value'] * LumB.rs[g,'cancer.p.value'] * Her2.rs[g,'cancer.p.value'] * Basal.rs[g,'cancer.p.value'])
# #   p.value <- 1- pchisq(s,df = 2 * 4)
# # }
# # names(combined.p.value.vec) <- marker.gene.symbol
# # combined.p.value.vec        <- sort(combined.p.value.vec)
####################### compute immune cell infiltration score ###########################

# flag                                  <- liver.vs.breast.immune.cell.signature$ensemble.id %in% rownames(cancer.data)
# liver.vs.breast.immune.cell.signature <- liver.vs.breast.immune.cell.signature[flag,]
# cell.type.vec                         <- unique(liver.vs.breast.immune.cell.signature$cell.type)
# immune.score.matrix <- foreach(cell.type = cell.type.vec,.combine='rbind') %do% {
#     g <- liver.vs.breast.immune.cell.signature$ensemble.id[liver.vs.breast.immune.cell.signature$cell.type == cell.type] %>% as.character  
#     if(length(g) > 1){
#         x <- apply(cancer.data[g,],2,function(x) 2^mean(x))
#         c(x)
#     }
#     else{
#         2^cancer.data[g,]
#     }
# }
# rownames(immune.score.matrix) <- cell.type.vec



