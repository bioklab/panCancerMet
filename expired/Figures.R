


require(ggplot2)
source('client-side/code/Manuscript/ggplot.style.R')
load('client-side/output/DE.R.output/DE.RData')

g  <- intersect(rownames(de.res.liver.vs.breast),rownames(de.res.metastasis.liver.vs.breast.lumb))
df <- data.frame(y=de.res.liver.vs.breast[g,'log2FoldChange'],x=de.res.metastasis.liver.vs.breast.lumb[g,'log2FoldChange'])
rownames(df) <- g
ggplot(df,aes(x=x,y=y)) + geom_point() + ggplot.style + xlab('M.vs.P') + ylab('Liver.vs.Breast') + xlim(c(-15,15)) + ylim(c(-15,15)) + geom_abline(intercept = 0,slope=1)




####
gene <- 'ENSG00000137100' # outlier in LumB



res <- de.res.liver.vs.breast
de.res.liver.vs.breast.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.01]
de.res.liver.vs.breast.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.01]



res <- de.res.metastasis.liver.vs.breast.basal
de.res.metastasis.liver.vs.breast.basal.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.01]
de.res.metastasis.liver.vs.breast.basal.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.01]



res <- de.res.metastasis.liver.vs.breast.her2
de.res.metastasis.liver.vs.breast.her2.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.01]
de.res.metastasis.liver.vs.breast.her2.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.01]

res <- de.res.metastasis.liver.vs.breast.luma
de.res.metastasis.liver.vs.breast.luma.up.gene <- rownames(res)[res$log2FoldChange > 1 & res$padj < 0.01]
de.res.metastasis.liver.vs.breast.luma.dn.gene <- rownames(res)[res$log2FoldChange < -1 & res$padj < 0.01]



g <- intersect(rownames(de.res.liver.vs.breast),rownames(de.res.metastasis.liver.vs.breast.her2))
plot(x=de.res.liver.vs.breast[g,2],y=de.res.metastasis.liver.vs.breast.her2[g,2],xlim=c(-15,15),ylim=c(-15,15),xlab='liver.vs.breast',ylab='metastatic.vs.primary')
abline(v=c(-1.5,1.5))
abline(h=c(-1.5,1.5))
lines(c(-15,15),c(-15,15))
dd <- data.frame(x=de.res.liver.vs.breast[g,2],y=de.res.metastasis.liver.vs.breast.her2[g,2])
rownames(dd) <- g
View(dd[dd$x > 1.5 & dd$y < -1.5,])
View(dd[dd$x < -1.5 & dd$y > 1.5,])



g <- intersect(rownames(de.res.liver.vs.breast),rownames(de.res.metastasis.liver.vs.breast.basal))
plot(x=de.res.liver.vs.breast[g,2],y=de.res.metastasis.liver.vs.breast.basal[g,2],xlim=c(-15,15),ylim=c(-10,10))
abline(v=c(-1.5,1.5))
abline(h=c(-1.5,1.5))
dd <- data.frame(x=de.res.liver.vs.breast[g,2],y=de.res.metastasis.liver.vs.breast.basal[g,2])
rownames(dd) <- g
View(dd[dd$x > 1.5 & dd$y < -1.5,])
View(dd[dd$x < -1.5 & dd$y > 1.5,])



MET500.sample <- MET500.breast.cancer.polyA.LumB.sample
MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
TCGA.sample   <- TCGA.breast.cancer.polyA.LumB.sample

rs.MET500 <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
rs.TCGA <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)

#View(rs.MET500$correlation.matrix)
s1 <-  rs.MET500$correlation.matrix[,'EFM192A_BREAST']

s2 <-  rs.TCGA$correlation.matrix[,'EFM192A_BREAST']

# plot(x=s1,y=MET500.log2.fpkm.matrix['ENSG00000153563',MET500.sample])
# 
# GATA2 <- 'ENSG00000179348'
# GATAD2A <- 'ENSG00000167491'
y <- c(MET500.log2.fpkm.matrix['ENSG00000136826',MET500.sample],TCGA.breast.cancer.log2.fpkm.matrix['ENSG00000136826',TCGA.sample])
x <- c(s1,s2)
plot(x,y)

color.vec <- ifelse(names(y) %in% names(s1),'red','black')
plot(x,y,col=color.vec,pch=19)

ADAT3 <-  'ENSG00000213638'


##############
shuffle.gene.set <- foreach(x=xCell.corrected.gene.set) %do% {
  sample(rownames(cancer.data.for.xCell),length(x)) 
  
}


shuffle.cancer.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],
                                           shuffle.gene.set, method = "ssgsea",
                                   ssgsea.norm = FALSE)
x             <- c(purity.TCGA,purity.MET500)

shuffle.cancer.de.df <- foreach(i =1:nrow(cancer.ssgsea.scores),.combine='rbind') %do% {
  df                <- data.frame(x=x,y=shuffle.cancer.ssgsea.scores[i,])
  loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
  adjusted.expr     <- loess.fit$residuals
  p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
  delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
  data.frame(cancer.delta=delta,cancer.p.value=p.value)
  
}
rownames(shuffle.cancer.de.df)  <- rownames(cancer.ssgsea.scores)
#cancer.de.df$cancer.fdr <- p.adjust(cancer.de.df$cancer.p.value,method='fdr')



xCell.data$signatures[[1]]@geneIds




data                  <- BRACA.log2.fpkm.matrix
data                  <- data[rownames(data) %in% mapping.df$ENSEMBL,]
rownames(data)        <- mapping.df[rownames(data),'SYMBOL']
hehe <- GSVA::gsva(expr=data,
                   xCell.corrected.gene.set, method = "ssgsea",
                   ssgsea.norm = FALSE)





MET500.sample <- MET500.breast.cancer.polyA.Basal.sample
MET500.sample <- intersect(MET500.sample,MET500.liver.sample)
TCGA.sample   <- TCGA.breast.cancer.polyA.Basal.sample

cancer.ssgsea.scores <- GSVA::gsva(expr=cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],
                                   xCell.data$signatures, method = "ssgsea",
                                   ssgsea.norm = FALSE) # hmm, maybe here should set some cutoff, ssgsea score smaller than this cut-off truncate to zero. Let us check CCLE data

#xCell.score.matrix <- xCellAnalysis(cancer.data.for.xCell[,c(TCGA.sample,MET500.sample)],file.name = 'client-side/output/TME.R.output/basal.raw.score.txt' )

rs.MET500     <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.sample],          expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
rs.TCGA       <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
purity.MET500 <- rs.MET500$correlation.matrix[,'HCC70_BREAST'] # use HCC70 cell line
purity.TCGA   <- rs.TCGA$correlation.matrix[,'HCC70_BREAST']
x             <- c(purity.TCGA,purity.MET500)


cancer.geneset.da.df <- foreach(i =1:nrow(cancer.ssgsea.scores),.combine='rbind') %do% {
  df                <- data.frame(x=x,y=cancer.ssgsea.scores[i,])
  loess.fit         <- loess(data = df,formula= y ~ x,span=0.3,degree=1,family="symmetric",iterations=4,surface="direct") # see https://stat.ethz.ch/pipermail/bioconductor/2003-September/002337.html
  adjusted.expr     <- loess.fit$residuals
  p.value           <- wilcox.test(adjusted.expr[MET500.sample],adjusted.expr[TCGA.sample])$p.value
  delta             <- median(adjusted.expr[MET500.sample]) -  median(adjusted.expr[TCGA.sample])
  data.frame(cancer.delta=delta,cancer.p.value=p.value)
  
}
rownames(cancer.geneset.da.df)  <- rownames(cancer.ssgsea.scores)
cancer.geneset.da.df$cancer.fdr <- p.adjust(cancer.geneset.da.df$cancer.p.value,method='fdr')
da.gene.set                     <- rownames(cancer.geneset.da.df)[cancer.geneset.da.df$cancer.fdr < 0.05]

Basal.cancer.geneset.da.df      <- cancer.geneset.da.df
Basal.da.gene.set               <- da.gene.set








