# associaiton between staging and DE genes
# survival analysis for DE genes
load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
load('server-side/RData//Breast Invasive Carcinoma.RData')
load('~/Project/InSilicoCRISPR/client-side/output/organize.TCGA.clinical.data.R.output/organize.TCGA.clinical.data.RData')
load('client-side/output/tumor.purity.based.on.cell.line.R.output//tumor.purity.based.on.cell.line.RData')
source('client-side/code/util.R')
require(ggplot2)
require(dplyr)
require(GSVA)
require(survival)


######## DE genes, survival analysis ###################

TCGA.breast.cancer.log2.read.count.matrix <- log2.read.count.matrix
TCGA.breast.cancer.log2.fpkm.matrix       <- log2.fpkm.matrix
protein.coding.gene.id                    <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character


perform.survival.analysis <- function(method= 'continuous') {
    adjusted.expr.matrix  <- TCGA.breast.cancer.log2.fpkm.matrix[target.gene,subtype.sample]
    clinical.data         <- TCGA.clinical.data.list$BRCA

    rs.df <- foreach(g=rownames(adjusted.expr.matrix),.combine='rbind') %do% {
        adjusted.expr  <- adjusted.expr.matrix[g,]
        patient.name   <- gsub(adjusted.expr %>% names,pattern = '\\-01',replacement = '')
        df <- data.frame(expr=adjusted.expr,
                         status=clinical.data[patient.name,'vital_status'], 
                         d1=clinical.data[patient.name,'days_to_death'] ,
                         d2=clinical.data[patient.name,'days_to_last_followup'],
                         age= clinical.data[patient.name,'years_to_birth'],
                         pathologic_stage= clinical.data[patient.name,'pathologic_stage']
                         )
        df[is.na(df$d1),'d1'] <- 0
        df[is.na(df$d2),'d2'] <- 0
        df$time <- (df$d1 + df$d2) / 30
        df      <- df[complete.cases(df),]
        if(method == 'group'){
            q                          <- quantile(df$expr)
            df$group                   <- 'median'
            df[df$expr > q[3],'group'] <- 'high'
            df[df$expr <= q[3],'group']<- 'low'
            df                         <- df[df$group!= 'median',]
            df$group                   <- factor(df$group,levels = c('low','high'))
            sdf                        <- survdiff(Surv(time, status) ~ group , data=df)
            p.val                      <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
            hr                         <- (sdf$obs[2]/sdf$exp[2]) / (sdf$obs[1]/sdf$exp[1])
            data.frame(p.value=p.val,effect.size=hr %>% log2)
        }else{
            cox                        <- coxph(Surv(time, status) ~ expr + age  , data = df)
            x                          <- summary(cox)
            data.frame(p.value=x$coefficients['expr','Pr(>|z|)'],effect.size=x$coefficients['expr','coef'])
        }
    }
    rownames(rs.df) <- rownames(adjusted.expr.matrix)
    rs.df
}


subtype.sample                   <- intersect(pure.TCGA.breast.cancer.polyA.Basal.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
target.gene                      <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
basal.survival.rs.up             <- perform.survival.analysis(method='group')
target.gene                      <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
basal.survival.rs.dn             <- perform.survival.analysis(method='group')


subtype.sample        <- intersect(pure.TCGA.breast.cancer.polyA.LumB.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
target.gene           <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.up.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
lumb.survival.rs.up   <- perform.survival.analysis(method='group')
target.gene           <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
lumb.survival.rs.dn   <- perform.survival.analysis(method='group')


subtype.sample      <- intersect(pure.TCGA.breast.cancer.polyA.Her2.sample,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
target.gene         <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output//her2.up.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
her2.survival.rs.up <- perform.survival.analysis(method='group')
target.gene         <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
her2.survival.rs.dn <- perform.survival.analysis(method='group')








############### DE genes - stage  assocaition analysis #########################
clinical.data         <- TCGA.clinical.data.list$BRCA
stage.vec             <- clinical.data$pathologic_stage %>% as.character
names(stage.vec)      <- paste(rownames(clinical.data),'-01',sep='')
flag                  <- stage.vec %in% c('stage i','stage ia','stage ib')
stage.vec[flag]       <- 'stage I'
flag                  <- stage.vec %in% c('stage ii','stage iia','stage iib')
stage.vec[flag]       <- 'stage II'
flag                  <- stage.vec %in% c('stage iiic','stage iiia','stage iiib','stage iii')
stage.vec[flag]       <- 'stage III'
flag                  <- stage.vec %in% c('stage iv')
stage.vec[flag]       <- 'stage IV'
flag                  <- stage.vec %in% c('stage x')
stage.vec[flag]       <- 'stage X'
stage.vec             <- stage.vec[stage.vec != 'stage X']


age.vec               <- clinical.data$years_to_birth
names(age.vec)        <- paste(rownames(clinical.data),'-01',sep='')


require(GSVA)
perform.stage.expression.analysis <- function(){
    g                   <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,subtype.sample])
    ssgsea.score.matrix <- gsva(expr = TCGA.breast.cancer.log2.fpkm.matrix[g,subtype.sample],
                                method='ssgsea',
                                ssgsea.norm=FALSE,
                                gset.idx.list=list(up.gene = up.gene,dn.gene = dn.gene)
                                )
    gene.score.df <- data.frame(up.score= ssgsea.score.matrix['up.gene',],dn.score = ssgsea.score.matrix['dn.gene',],stage= stage.vec[subtype.sample])
    flag          <- gene.score.df$stage %in% c('stage I','stage II','stage III','stage IV')
    gene.score.df <- gene.score.df[flag,]
    #ggplot(gene.score.df,aes(x=stage,y=up.score)) + geom_boxplot()
    #ggplot(gene.score.df,aes(x=stage,y=dn.score)) + geom_boxplot()
    gene.score.df
    
}

########################  Basal subtype ########################
subtype.sample <- pure.TCGA.breast.cancer.polyA.Basal.sample
subtype.sample <- intersect(subtype.sample,names(stage.vec))
up.gene        <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
dn.gene        <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
basal.stage.association.rs <- perform.stage.expression.analysis()

########################  Her2 subtype ########################
subtype.sample <- pure.TCGA.breast.cancer.polyA.Her2.sample
subtype.sample <- intersect(subtype.sample,names(stage.vec))
up.gene  <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.up.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
dn.gene  <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
her2.stage.association.rs <- perform.stage.expression.analysis()


########################  LumB subtype ########################
subtype.sample <- pure.TCGA.breast.cancer.polyA.LumB.sample
subtype.sample <- intersect(subtype.sample,names(stage.vec))
up.gene  <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.up.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
dn.gene  <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.dn.gene.immune.excluded.csv", stringsAsFactors=FALSE)$x
lumb.stage.association.rs <- perform.stage.expression.analysis()

save(file = 'client-side/output/DE.gene.clinics.R.output/DE.gene.clinics.RData',list=c('lumb.stage.association.rs','her2.stage.association.rs','basal.stage.association.rs',
                                                                                       'lumb.survival.rs.up','lumb.survival.rs.dn',
                                                                                       'her2.survival.rs.up','her2.survival.rs.dn',
                                                                                       'basal.survival.rs.up','basal.survival.rs.dn'
                                                                                       )
)









############################# Trash code #################################################################



# load('~/Project/InSilicoCRISPR/client-side/output/organize.TCGA.mutation.data.R.output/organize.TCGA.mutation.data.RData')
# mutation.data <- TCGA.mutation.data.list$BRCA
# mutation.burden <- apply(mutation.data,1,sum)
# df$mutation.burden <- mutation.burden[rownames(df)]
# 
# ggplot(df)+geom_point(aes(x=purity,y=mutation.burden %>% log10,color=stage))
# ggplot(df)+geom_boxplot(aes(x=stage,y=mutation.burden %>% log,color=stage))
# 
# 
# 
# 
# 
# ##############
# df                 <- df[complete.cases(df),]
# df$mutation.burden <- log10(df$mutation.burden+1)
# h.sample           <- df$sample.id[df$mutation.burden > quantile(df$mutation.burden)[4]]
# l.sample           <- df$sample.id[df$mutation.burden < quantile(df$mutation.burden)[2]]
# g1                 <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,h.sample])
# g2                 <- get.expressed.gene(TCGA.breast.cancer.log2.fpkm.matrix[,l.sample])
# expressed.gene     <- intersect(c(g1,g2) %>% unique,protein.coding.gene.id)
# 
# 
# res.df <- foreach(g= expressed.gene,.combine='rbind') %do% {
#     effect.size <- median(TCGA.breast.cancer.log2.fpkm.matrix[g,h.sample]) -   median(TCGA.breast.cancer.log2.fpkm.matrix[g,l.sample])
#     p.value     <- wilcox.test(TCGA.breast.cancer.log2.fpkm.matrix[g,h.sample],TCGA.breast.cancer.log2.fpkm.matrix[g,l.sample])$p.value
#     data.frame(effect.size = effect.size,p.value=p.value)
# }
# rownames(res.df) <- expressed.gene
# plot(x=res.df$effect.size,y=-1 * (res.df$p.value %>% log10))
# res.df$fdr <- p.adjust(res.df$p.value,method='fdr')
# 
# up.gene <- rownames(res.df)[res.df$fdr < 0.1 & res.df$effect.size > 0.5]
# dn.gene <- rownames(res.df)[res.df$fdr < 0.1 & res.df$effect.size < -0.5]
# 
# up.gene.annotation <- AnnotationDbi::select(x = org.Hs.eg.db,keys = up.gene,keytype = 'ENSEMBL',columns=c('ENSEMBL','SYMBOL'))
# dn.gene.annotation <- AnnotationDbi::select(x = org.Hs.eg.db,keys = dn.gene,keytype = 'ENSEMBL',columns=c('ENSEMBL','SYMBOL'))
# 
# 
# 
# 
# df$hyper.mutation <- 'NULL'
# df[h.sample,'hyper.mutation'] <- 'high'
# df[l.sample,'hyper.mutation'] <- 'low'
# 
# 
# s              <- c(h.sample,l.sample[1:26])
# expr.matrix <- TCGA.breast.cancer.log2.read.count.matrix[expressed.gene,s]
# expr.matrix <- 2^expr.matrix - 1
# 
# 
# dds           <- DESeqDataSetFromMatrix(countData = round(expr.matrix),
#                                         colData = df[s,],
#                                         design = ~ hyper.mutation + purity)
# dds        <- DESeq(dds)
# 
# res <- results(dds,contrast = c('hyper.mutation','high','low'))   %>% as.data.frame
# res <- res[complete.cases(res),]
# res <- res[order(res$padj),]
# 
# plot(x=res$log2FoldChange,y=-1 * (res$pvalue %>% log10))
# 
# View(res)
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
# 
# 
# 
# 
# 
# 
# 
# 
# 
# plot(x=res.2.vs.1[up.gene,'log2FoldChange'],y=res.3.vs.2[up.gene,'log2FoldChange'])
# plot(x=res.2.vs.1[dn.gene,'log2FoldChange'],y=res.3.vs.2[dn.gene,'log2FoldChange'])
# 
# HES4 <- 'ENSG00000188290'
# 
# 
# 
# 
# 
# Basal.up.gene          <- read.csv(file='client-side/output/DE.breast.cancer.R.output/basal.up.csv')$x %>% as.character
# Basal.dn.gene          <- read.csv(file='client-side/output/DE.breast.cancer.R.output/basal.dn.csv')$x %>% as.character
# Basal.log2.fpkm.matrix <- TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.Basal.sample]
# 
# 
# 
# 
# score.matrix <- gsva(expr=Basal.log2.fpkm.matrix,gset.idx.list = list(up.gene=Basal.up.gene,dn.gene=Basal.dn.gene),method='ssgsea',ssgsea.norm=FALSE)
# 
# df <- data.frame(score=score.matrix %>% t,stage=stage.vec[colnames(score.matrix)],purity=tumor.purity.based.on.cell.line.vec[colnames(score.matrix)])
# df <- df[df$stage %in% c('stage I','stage II','stage III'),]
# 
# ggplot(df) + geom_boxplot(aes(x=stage,y=score.up.gene)) + geom_point(aes(x=stage,y=score.up.gene)) + ggplot.style
# ggplot(df) + geom_boxplot(aes(x=stage,y=score.dn.gene)) + geom_point(aes(x=stage,y=score.dn.gene)) + ggplot.style
# 
# wilcox.test(df$score.up.gene[df$stage == 'stage I'],df$score.up.gene[df$stage == 'stage II'])
# wilcox.test(df$score.up.gene[df$stage == 'stage II'],df$score.up.gene[df$stage == 'stage III'])
# wilcox.test(df$score.up.gene[df$stage == 'stage I'],df$score.up.gene[df$stage == 'stage III'])
# wilcox.test(df$score.dn.gene[df$stage == 'stage I'],df$score.dn.gene[df$stage == 'stage II'])
# wilcox.test(df$score.dn.gene[df$stage == 'stage II'],df$score.dn.gene[df$stage == 'stage III'])
# wilcox.test(df$score.dn.gene[df$stage == 'stage I'],df$score.dn.gene[df$stage == 'stage III'])
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
# ########################  LumB subtype ########################
# LumB.up.gene          <- read.csv(file='client-side/output/DE.breast.cancer.R.output/lumb.up.csv')$x %>% as.character
# LumB.dn.gene          <- read.csv(file='client-side/output/DE.breast.cancer.R.output/lumb.dn.csv')$x %>% as.character
# LumB.log2.fpkm.matrix <- TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.LumB.sample]
# 
# clinical.data         <- TCGA.clinical.data.list$BRCA
# stage.vec             <- clinical.data$pathologic_stage %>% as.character
# names(stage.vec)      <- paste(rownames(clinical.data),'-01',sep='')
# stage.vec             <- stage.vec[names(stage.vec) %in% colnames(TCGA.breast.cancer.log2.fpkm.matrix)]
# stage.vec             <- stage.vec[names(stage.vec) %in%  TCGA.breast.cancer.polyA.LumB.sample]
# 
# flag <- stage.vec %in% c('stage i','stage ia')
# stage.vec[flag] <- 'stage I'
# flag <- stage.vec %in% c('stage ii','stage iia','stage iib')
# stage.vec[flag] <- 'stage II'
# flag <- stage.vec %in% c('stage iiic','stage iiia','stage iiib')
# stage.vec[flag] <- 'stage III'
# 
# 
# score.matrix <- gsva(expr=LumB.log2.fpkm.matrix,gset.idx.list = list(up.gene=LumB.up.gene,dn.gene=LumB.dn.gene),method='ssgsea',ssgsea.norm=FALSE)
# 
# df <- data.frame(score=score.matrix %>% t,stage=stage.vec[colnames(score.matrix)],purity=tumor.purity.based.on.cell.line.vec[colnames(score.matrix)])
# df <- df[df$stage %in% c('stage I','stage II','stage III'),]
# 
# ggplot(df) + geom_boxplot(aes(x=stage,y=score.up.gene)) + geom_point(aes(x=stage,y=score.up.gene)) + ggplot.style
# ggplot(df) + geom_boxplot(aes(x=stage,y=score.dn.gene)) + geom_point(aes(x=stage,y=score.dn.gene)) + ggplot.style
# 
# 
# ########################  LumA subtype ########################
# LumA.up.gene          <- read.csv(file='client-side/output/DE.breast.cancer.R.output/luma.up.csv')$x %>% as.character
# LumA.dn.gene          <- read.csv(file='client-side/output/DE.breast.cancer.R.output/luma.dn.csv')$x %>% as.character
# LumA.log2.fpkm.matrix <- TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.LumA.sample]
# 
# clinical.data         <- TCGA.clinical.data.list$BRCA
# stage.vec             <- clinical.data$pathologic_stage %>% as.character
# names(stage.vec)      <- paste(rownames(clinical.data),'-01',sep='')
# stage.vec             <- stage.vec[names(stage.vec) %in% colnames(TCGA.breast.cancer.log2.fpkm.matrix)]
# stage.vec             <- stage.vec[names(stage.vec) %in%  TCGA.breast.cancer.polyA.LumA.sample]
# 
# flag <- stage.vec %in% c('stage i','stage ia')
# stage.vec[flag] <- 'stage I'
# flag <- stage.vec %in% c('stage ii','stage iia','stage iib')
# stage.vec[flag] <- 'stage II'
# flag <- stage.vec %in% c('stage iiic','stage iiia','stage iiib')
# stage.vec[flag] <- 'stage III'
# 
# 
# score.matrix <- gsva(expr=LumA.log2.fpkm.matrix,gset.idx.list = list(up.gene=LumA.up.gene,dn.gene=LumA.dn.gene),method='ssgsea',ssgsea.norm=FALSE)
# 
# df <- data.frame(score=score.matrix %>% t,stage=stage.vec[colnames(score.matrix)],purity=tumor.purity.based.on.cell.line.vec[colnames(score.matrix)])
# df <- df[df$stage %in% c('stage I','stage II','stage III'),]
# 
# ggplot(df) + geom_boxplot(aes(x=stage,y=score.up.gene)) + geom_point(aes(x=stage,y=score.up.gene)) + ggplot.style
# ggplot(df) + geom_boxplot(aes(x=stage,y=score.dn.gene)) + geom_point(aes(x=stage,y=score.dn.gene)) + ggplot.style
# 
# 
# wilcox.test(df$score.up.gene[df$stage == 'stage I'],df$score.up.gene[df$stage == 'stage II'])
# wilcox.test(df$score.up.gene[df$stage == 'stage II'],df$score.up.gene[df$stage == 'stage III'])
# wilcox.test(df$score.up.gene[df$stage == 'stage I'],df$score.up.gene[df$stage == 'stage III'])
# wilcox.test(df$score.dn.gene[df$stage == 'stage I'],df$score.dn.gene[df$stage == 'stage II'])
# wilcox.test(df$score.dn.gene[df$stage == 'stage II'],df$score.dn.gene[df$stage == 'stage III'])
# wilcox.test(df$score.dn.gene[df$stage == 'stage I'],df$score.dn.gene[df$stage == 'stage III'])
# 
# 
# 
# 
# ########################  Her2 subtype ########################
# Her2.up.gene          <- read.csv(file='client-side/output/DE.breast.cancer.R.output/her2.up.csv')$x %>% as.character
# Her2.dn.gene          <- read.csv(file='client-side/output/DE.breast.cancer.R.output/her2.dn.csv')$x %>% as.character
# Her2.log2.fpkm.matrix <- TCGA.breast.cancer.log2.fpkm.matrix[,TCGA.breast.cancer.polyA.Her2.sample]
# 
# clinical.data         <- TCGA.clinical.data.list$BRCA
# stage.vec             <- clinical.data$pathologic_stage %>% as.character
# names(stage.vec)      <- paste(rownames(clinical.data),'-01',sep='')
# stage.vec             <- stage.vec[names(stage.vec) %in% colnames(TCGA.breast.cancer.log2.fpkm.matrix)]
# stage.vec             <- stage.vec[names(stage.vec) %in%  TCGA.breast.cancer.polyA.Her2.sample]
# 
# flag <- stage.vec %in% c('stage i','stage ia')
# stage.vec[flag] <- 'stage I'
# flag <- stage.vec %in% c('stage ii','stage iia','stage iib')
# stage.vec[flag] <- 'stage II'
# flag <- stage.vec %in% c('stage iiic','stage iiia','stage iiib')
# stage.vec[flag] <- 'stage III'
# 
# 
# score.matrix <- gsva(expr=Her2.log2.fpkm.matrix,gset.idx.list = list(up.gene=Her2.up.gene,dn.gene=Her2.dn.gene),method='ssgsea',ssgsea.norm=FALSE)
# 
# df <- data.frame(score=score.matrix %>% t,stage=stage.vec[colnames(score.matrix)],purity=tumor.purity.based.on.cell.line.vec[colnames(score.matrix)])
# df <- df[df$stage %in% c('stage I','stage II','stage III'),]
# 
# ggplot(df) + geom_boxplot(aes(x=stage,y=score.up.gene)) + geom_point(aes(x=stage,y=score.up.gene)) + ggplot.style
# ggplot(df) + geom_boxplot(aes(x=stage,y=score.dn.gene)) + geom_point(aes(x=stage,y=score.dn.gene)) + ggplot.style
# 
# 
# wilcox.test(df$score.up.gene[df$stage == 'stage I'],df$score.up.gene[df$stage == 'stage II'])
# wilcox.test(df$score.up.gene[df$stage == 'stage II'],df$score.up.gene[df$stage == 'stage III'])
# wilcox.test(df$score.up.gene[df$stage == 'stage I'],df$score.up.gene[df$stage == 'stage III'])
# wilcox.test(df$score.dn.gene[df$stage == 'stage I'],df$score.dn.gene[df$stage == 'stage II'])
# wilcox.test(df$score.dn.gene[df$stage == 'stage II'],df$score.dn.gene[df$stage == 'stage III'])
# wilcox.test(df$score.dn.gene[df$stage == 'stage I'],df$score.dn.gene[df$stage == 'stage III'])
# 






# NOTCH1  <- 'ENSG00000148400'
# NOTCH3  <- 'ENSG00000074181'
# HES4    <- 'ENSG00000188290' # wow, bimodal distribution in Basal
# PTPN2   <- 'ENSG00000175354'
# KLF4    <- 'ENSG00000136826'
# MYC     <- 'ENSG00000136997'
# GATAD2A <- 'ENSG00000167491'
# BRD4    <- 'ENSG00000141867'
# TSSK6   <- 'ENSG00000178093' #hmm, seems interesting
# 
# g <- HES4
# 
# df <- data.frame(stage=stage.vec,expr = TCGA.breast.cancer.log2.fpkm.matrix[g,names(stage.vec)])
# ggplot(df) + geom_boxplot(aes(x=stage,y=expr)) + geom_point(aes(x=stage,y=expr))
# wilcox.test(df$expr[df$stage == 'stage I'],df$expr[df$stage == 'stage II'])
# 
# load("/Users/liuke/Project/InSilicoCRISPR/client-side/output/compute.TCGA.CD8.T.cell.level.R.output/compute.TCGA.CD8.T.cell.level.RData")
# 
# df        <- df[rownames(df) %in% rownames(TCGA.CD8.T.cell.level.df),]
# df        <- df[df$stage %in% c('stage iia','stage iib'),]
# df$purity <- TCGA.CD8.T.cell.level.df[rownames(df),'tumor.purity']
# ggplot(df,aes(x=purity,y=expr)) + geom_point() + stat_smooth(method='lm')
# 
# 
# 
# 
# 
# 
# 
# load('~/Project/InSilicoCRISPR/client-side/output/organize.TCGA.mutation.data.R.output/organize.TCGA.mutation.data.RData')
# BRCA.mutation   <- TCGA.mutation.data.list$BRCA
# mutation.burden <- apply(BRCA.mutation,1,sum)
# mutation.burden <- mutation.burden[names(mutation.burden) %in% colnames(TCGA.breast.cancer.log2.fpkm.matrix)]
# 
# # age.vec                                <- clinical.data$years_to_birth 
# # names(age.vec)                         <- paste(rownames(clinical.data),'-01',sep='')
# 
# 
# df <- data.frame(mutation.burden=mutation.burden,expr = TCGA.breast.cancer.log2.fpkm.matrix[NOTCH3,names(mutation.burden)])
# ggplot(df,aes(x=log10(mutation.burden) ,y=expr)) + geom_point()
