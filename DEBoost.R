require(DESeq2)
require(dplyr)
require(foreach)
require(segmented)
load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')


########## get expressed genes ##############
get.expressed.gene <- function(expr.matrix,cut.off=1){
    m.expr <- apply(expr.matrix,1,median)
    rownames(expr.matrix)[m.expr >= cut.off ]
}


########## DE analysis ##############
perform.DE.analysis.between.TRE.and.CON <- function(CON.log2.tpm.matrix, CON.log2.read.count.matrix, TRE.log2.tpm.matrix,TRE.log2.read.count.matrix)
{
    g1             <- get.expressed.gene(CON.log2.tpm.matrix)
    g2             <- get.expressed.gene(TRE.log2.tpm.matrix)
    expressed.gene <- c(g1,g2) %>% unique
  
    read.count.matrix    <- cbind(TRE.log2.read.count.matrix[expressed.gene,],CON.log2.read.count.matrix[expressed.gene,])
    read.count.matrix    <- 2^read.count.matrix - 1
  
    df             <- data.frame(condition=c(rep(x='TRE',times=TRE.log2.tpm.matrix %>% ncol), rep(x='CON',times=CON.log2.tpm.matrix %>% ncol) ))
    df$condition   <- factor(df$condition,levels = c('CON','TRE'))
  
    dds           <- DESeqDataSetFromMatrix(countData = round(read.count.matrix),
                                            colData = df,
                                            design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds,contrast = c('condition','TRE','CON'),cooksCutoff = FALSE) %>% as.data.frame
    res$padj <- p.adjust(res$pvalue,method='fdr')
    #res <- res[complete.cases(res),]     
    res <- res[order(res$pvalue),]
    res
}

########## DEBoost filtering ##############

DEBoost.filtering <- function(deseq2.M.vs.P.res, deseq2.R.vs.P.res, MET.log2.tpm.matrix, PRI.log2.tpm.matrix, REF.log2.tpm.matrix, TCGA.best.cell.line,MET500.best.cell.line ){
    up.gene            <- rownames(deseq2.M.vs.P.res)[deseq2.M.vs.P.res$log2FoldChange > 1  & deseq2.M.vs.P.res$padj < 0.05]
    dn.gene            <- rownames(deseq2.M.vs.P.res)[deseq2.M.vs.P.res$log2FoldChange < -1 & deseq2.M.vs.P.res$padj < 0.05]
  
    ######### piecewise linear regression to identify ref-tissue specific genes #########
    psi     <- -1
    while(psi < 0){
        c.gene        <- intersect(rownames(deseq2.M.vs.P.res),rownames(deseq2.R.vs.P.res))
        x             <- deseq2.R.vs.P.res[c.gene,'log2FoldChange']
        y             <- deseq2.M.vs.P.res[c.gene,'log2FoldChange']
        lin.mod       <- lm(y~x)
        segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)    
        tmp           <- summary(segmented.mod)
        psi           <- tmp$psi[1,'Est.']
    }
    
    ref.specific.gene  <- rownames(deseq2.R.vs.P.res)[deseq2.R.vs.P.res$log2FoldChange > psi]
    
    r           <- segmented.mod$residuals
    names(r)    <- c.gene
    flag        <- x > psi
    r           <- r[flag]
    r.scale     <- scale(r) 
    up.outlier.gene <- rownames(r.scale)[r.scale[,1] > 3 ]
    dn.outlier.gene <- rownames(r.scale)[r.scale[,1] < -3 ]
    liver.marker.gene               <- ref.specific.gene
    MET.marker.gene.expr.matrix     <- 2^ MET.log2.tpm.matrix[liver.marker.gene,] - 1
    REF.marker.gene.expr            <- apply(2^REF.log2.tpm.matrix[liver.marker.gene,] - 1, 1, median)
    alpha.vec                       <- foreach(i=1:ncol(MET.marker.gene.expr.matrix),.combine='c') %do% {
        x <- log2(MET.marker.gene.expr.matrix[,i] / REF.marker.gene.expr) %>% median 
        x
    }
    names(alpha.vec) <- colnames(MET.marker.gene.expr.matrix)
    
    most.pure.sample <- names(alpha.vec) [alpha.vec < log2(0.01) ]
    if(length(most.pure.sample) >= 2 & length(up.outlier.gene) >= 2) {
        up.outlier.gene.expr.0.01 <- apply(MET.log2.tpm.matrix[up.outlier.gene,most.pure.sample],1,median)
    }else{
        up.outlier.gene.expr.0.01 <- c(-1)
    }
    
    most.pure.sample <- names(alpha.vec) [alpha.vec < log2(0.05) ]
    if(length(most.pure.sample) >= 2 & length(up.outlier.gene) >= 2) {
        up.outlier.gene.expr.0.05 <- apply(MET.log2.tpm.matrix[up.outlier.gene,most.pure.sample],1,median)
    }else{
        up.outlier.gene.expr.0.05 <- c(-1)
    }
  
    # ####### gene filtering based on expression level #############
    # immune.gene.list        <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/immune.gene.list.csv", stringsAsFactors=FALSE)$x
    # PRI.median              <- apply(PRI.log2.tpm.matrix, 1, median)
    # MET.median              <- apply(MET.log2.tpm.matrix, 1, median)
    # median.max              <- apply(rbind(PRI.median,MET.median),2,max)
    # x                       <- median.max[immune.gene.list]
    # x                       <- x[is.na(x) == FALSE]
    # q                       <- quantile(x)
    # tpm.cut.off             <- q['75%'] + 1.5 * IQR(x)
    # #tpm.cut.off            <- log2(5 + 1)
    # 
    # low.expr.gene           <- names(median.max)[median.max < tpm.cut.off]
  
  
    cell.line.expr.max       <- apply(CCLE.log2.rpkm.matrix[,c(TCGA.best.cell.line,MET500.best.cell.line)],1,max)
    cell.line.expressed.gene <- names(cell.line.expr.max)[cell.line.expr.max >= log2(5 + 1)]
    
    
    
    ####### gene filtering based on wilcoxon rank test #######################
    wilcox.test.p.value.vec <- foreach(g=c(up.gene,dn.gene),.combine='c') %do% {
        wilcox.test(MET.log2.tpm.matrix[g,],PRI.log2.tpm.matrix[g,])$p.value  
    }
    names(wilcox.test.p.value.vec) <- c(up.gene,dn.gene)
    not.robust.gene                <- names(wilcox.test.p.value.vec) [ wilcox.test.p.value.vec > 0.05]
  
  
    uDE.gene <- c(ref.specific.gene,not.robust.gene) %>% unique()
    up.gene  <- setdiff(up.gene,uDE.gene)
    dn.gene  <- setdiff(dn.gene,uDE.gene)
    up.gene  <- intersect(up.gene,cell.line.expressed.gene)
    dn.gene  <- intersect(dn.gene,cell.line.expressed.gene)
    
  
    list(up.gene=c(up.gene),dn.gene=c(dn.gene),up.outlier.gene = up.outlier.gene,dn.outlier.gene = dn.outlier.gene,up.outlier.gene.expr.0.01 = up.outlier.gene.expr.0.01,up.outlier.gene.expr.0.05 = up.outlier.gene.expr.0.05)
  
  ######### gene filtering based on correlation with purity#########
  # p              <- tumor.purity.based.on.cell.line.vec[TCGA.sample]
  # cor.vec        <- cor(TCGA.breast.cancer.log2.tpm.matrix[,TCGA.sample] %>% t,p,method='spearman') %>% c
  # names(cor.vec) <- rownames(TCGA.breast.cancer.log2.tpm.matrix)
  # cor.data       <- cor.vec[c(up.gene,dn.gene)]
  # model          <- normalmixEM(cor.data)
  # mu.vec         <- model$mu
  # sigma.vec      <- model$sigma
  # if(mu.vec[1] < mu.vec[2]){
  #     mu     <- mu.vec[1]
  #     sigma  <- sigma.vec[1]
  # } else {
  #     mu      <- mu.vec[2]
  #     sigma   <- sigma.vec[2]
  # }
  # cor.cut.off  <- mu + 3.0 * sigma
  # neg.cor.gene <- names(cor.data)[ cor.data <= cor.cut.off]
  # ggplot.graph <- ggplot(data.frame(x=cor.data), aes(x=x)) + geom_histogram() + geom_density(aes(y=..density.. * 50)) + geom_vline(xintercept = cor.cut.off)
  # 
}



perform.DE.analysis.between.primary.and.metastatic.cancer <- function(PRI.log2.tpm.matrix, PRI.log2.read.count.matrix,
                                                                      MET.log2.tpm.matrix, MET.log2.read.count.matrix,
                                                                      REF.log2.tpm.matrix, REF.log2.read.count.matrix,
                                                                      TCGA.best.cell.line, MET500.best.cell.line
                                                                      )

{
  deseq2.M.vs.P.res <- perform.DE.analysis.between.TRE.and.CON(CON.log2.tpm.matrix = PRI.log2.tpm.matrix, CON.log2.read.count.matrix = PRI.log2.read.count.matrix,
                                                               TRE.log2.tpm.matrix = MET.log2.tpm.matrix, TRE.log2.read.count.matrix = MET.log2.read.count.matrix
                                                               )
  
  deseq2.R.vs.P.res <- perform.DE.analysis.between.TRE.and.CON(CON.log2.tpm.matrix = PRI.log2.tpm.matrix, CON.log2.read.count.matrix = PRI.log2.read.count.matrix,
                                                               TRE.log2.tpm.matrix = REF.log2.tpm.matrix, TRE.log2.read.count.matrix = REF.log2.read.count.matrix
  )
  
  tumor.intrinsic.DE.gene.rs  <- DEBoost.filtering(deseq2.M.vs.P.res = deseq2.M.vs.P.res,deseq2.R.vs.P.res = deseq2.R.vs.P.res, 
                                                   MET.log2.tpm.matrix = MET.log2.tpm.matrix, PRI.log2.tpm.matrix = PRI.log2.tpm.matrix,
                                                   REF.log2.tpm.matrix = REF.log2.tpm.matrix,TCGA.best.cell.line,MET500.best.cell.line)
  
  list(deseq2.M.vs.P.res = deseq2.M.vs.P.res, deseq2.R.vs.P.res = deseq2.R.vs.P.res, tumor.intrinsic.DE.gene.rs= tumor.intrinsic.DE.gene.rs)
  
}


