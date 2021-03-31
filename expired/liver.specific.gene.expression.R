require(ggplot2)
require(dplyr)
require(segmented)
source('client-side/code/Manuscript/ggplot.style.R')
load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')


################## Luminal B subtype ###################
c.gene        <- intersect(rownames(de.res.liver.vs.breast.lumb),rownames(de.res.metastasis.liver.vs.breast.lumb))
x             <- de.res.liver.vs.breast.lumb[c.gene,'log2FoldChange']
y             <- de.res.metastasis.liver.vs.breast.lumb[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=de.res.liver.vs.breast.lumb[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.lumb[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
df$residual   <- df$y - df$fitted
lumb.liver.specific.gene.df    <- df[df$x > psi,]
lumb.liver.specific.gene.df$sr <- scale(lumb.liver.specific.gene.df$residual)
ggplot(lumb.liver.specific.gene.df,aes(x=x,y=sr)) + geom_point(size=2.5) + ggplot.style 
View(lumb.liver.specific.gene.df)



################## Her2 subtype ###################
c.gene        <- intersect(rownames(de.res.liver.vs.breast.her2),rownames(de.res.metastasis.liver.vs.breast.her2))
x             <- de.res.liver.vs.breast.her2[c.gene,'log2FoldChange']
y             <- de.res.metastasis.liver.vs.breast.her2[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=de.res.liver.vs.breast.her2[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.her2[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
df$residual   <- df$y - df$fitted
her2.liver.specific.gene.df    <- df[df$x > psi,]
her2.liver.specific.gene.df$sr <- scale(her2.liver.specific.gene.df$residual)
ggplot(her2.liver.specific.gene.df,aes(x=x,y=sr)) + geom_point(size=2.5) + ggplot.style 
View(her2.liver.specific.gene.df)



################## Basal subtype ###################
c.gene        <- intersect(rownames(de.res.liver.vs.breast.basal),rownames(de.res.metastasis.liver.vs.breast.basal))
x             <- de.res.liver.vs.breast.basal[c.gene,'log2FoldChange']
y             <- de.res.metastasis.liver.vs.breast.basal[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=de.res.liver.vs.breast.basal[c.gene,'log2FoldChange'],y=de.res.metastasis.liver.vs.breast.basal[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
df$residual   <- df$y - df$fitted
basal.liver.specific.gene.df    <- df[df$x > psi,]
basal.liver.specific.gene.df$sr <- scale(basal.liver.specific.gene.df$residual)
ggplot(basal.liver.specific.gene.df,aes(x=x,y=sr)) + geom_point(size=2.5) + ggplot.style 
View(basal.liver.specific.gene.df)





load('client-side/output/validation.of.confounding.R.output/validation.of.confounding.RData')

c.gene        <- intersect(rownames(SRP043470.de.res.liver.vs.breast),rownames(SRP043470.de.res.metastasis.liver.vs.breast))
x             <- SRP043470.de.res.liver.vs.breast[c.gene,'log2FoldChange']
y             <- SRP043470.de.res.metastasis.liver.vs.breast[c.gene,'log2FoldChange']
lin.mod       <- lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=2)
tmp           <- summary(segmented.mod)
psi           <- tmp$psi[1,'Est.']
df            <- data.frame(x=SRP043470.de.res.liver.vs.breast[c.gene,'log2FoldChange'],y=SRP043470.de.res.metastasis.liver.vs.breast[c.gene,'log2FoldChange'],fitted = fitted(segmented.mod))
rownames(df)  <- c.gene
df$residual   <- df$y - df$fitted
liver.specific.gene.df <- df[df$x > psi,]
liver.specific.gene.df$sr <- scale(liver.specific.gene.df$residual)

ggplot(liver.specific.gene.df,aes(x=x,y=sr)) + geom_point(size=2.5) + ggplot.style

