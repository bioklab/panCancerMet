require(dplyr)
require(stringr)

file.name <- "client-side/Data/JASPAR.txt"
conn      <- file(file.name,open="r")
line      <- readLines(conn)
get.TF.name <- function(x) {
  if(grepl(x = x,pattern='>')){
    l <- strsplit(x=x,split='\t')  %>% unlist
    toupper(l[2])
  }else{
    NA
  }
}

TF.list <- sapply(line,get.TF.name)
TF.list <- TF.list[is.na(TF.list) == FALSE]
names(TF.list) <- NULL
JASPAR.TF.list <- TF.list
close(conn)


file.name <- "client-side/Data/CISTROME.factor.line.txt"
conn      <- file(file.name,open="r")
line      <- readLines(conn)
tmp       <- str_extract_all(pattern="id=\"[:alnum:]+\"",string=line[1],simplify = TRUE)
get.TF.name <- function(x) {
  x <- str_remove_all(string = x,pattern="\"")
  x <- str_remove_all(string = x,pattern="id=")
  x
}
TF.list <- sapply(tmp,get.TF.name)
TF.list <- TF.list[is.na(TF.list) == FALSE]
names(TF.list) <- NULL
CISTROME.TF.list <- TF.list
close(conn)

TF.list <- c(CISTROME.TF.list,JASPAR.TF.list) %>% unique


# SNAI1 <- 'ENSG00000124216'
# SNAI2 <- 'ENSG00000019549'
# ZEB1  <- 'ENSG00000148516'
# ZEB2  <- 'ENSG00000169554'
# TWIST <- 'ENSG00000122691'
# CDH1  <- 'ENSG00000039068'
# EMT.TF <- c(SNAI1,SNAI2,ZEB1,ZEB2,TWIST,CDH1)

basal.up.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
basal.dn.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
intersect(basal.up.gene.symbol,TF.list)
intersect(basal.dn.gene.symbol,TF.list)

lumb.up.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.up.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
lumb.dn.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/lumb.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
intersect(lumb.up.gene.symbol,TF.list)
intersect(lumb.dn.gene.symbol,TF.list)

her2.up.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.up.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
her2.dn.gene.symbol <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/analyze.DE.gene.R.output/her2.dn.gene.immune.excluded.annotation.df.csv", stringsAsFactors=FALSE)$SYMBOL
intersect(her2.up.gene.symbol,TF.list)
intersect(her2.dn.gene.symbol,TF.list)


x1 <- intersect(basal.dn.gene.symbol,TF.list)
x2 <- intersect(lumb.dn.gene.symbol, TF.list)
x3 <- intersect(her2.dn.gene.symbol, TF.list)
tmp.dn.TF <- table(c(x1,x2,x3)) %>% as.data.frame



# x1 <- intersect(basal.up.gene.symbol,TF.list)
# x2 <- intersect(lumb.up.gene.symbol, TF.list)
# x3 <- intersect(her2.up.gene.symbol, TF.list)
# tmp.up.TF <- table(c(x1,x2,x3)) %>% as.data.frame


load('client-side/output/analyze.DE.gene.R.output/common.DE.gene.RData')
ECM.gene <- strsplit(x=c.dn.gene.BP.filtered$geneID[1],split = '/') %>% unlist

require(org.Hs.eg.db)
require(AnnotationDbi)
ECM.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'SYMBOL',columns = c('ENSEMBL','SYMBOL'),keys = ECM.gene) 
ECM.gene.df <- ECM.gene.df[ECM.gene.df$ENSEMBL != 'ENSG00000275365',]
ECM.gene.ensemble.id <- ECM.gene.df$ENSEMBL %>% as.character()

load('server-side/RData/Breast Invasive Carcinoma.RData')
load('client-side/output/TCGA.breast.cancer.meta.R.output/TCGA.breast.cancer.meta.RData')
PRRX1 <- 'ENSG00000116132'


other.gene <- setdiff(rownames(log2.fpkm.matrix),c(PRRX1,ECM.gene.ensemble.id))

Basal.PRRX1.ECM.cor.vec  <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.Basal.sample] ,log2.fpkm.matrix[ECM.gene.ensemble.id,pure.TCGA.breast.cancer.polyA.Basal.sample] %>% t,method='spearman') %>% c
Basal.PRRX1.other.cor.vec <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.Basal.sample] ,log2.fpkm.matrix[other.gene,pure.TCGA.breast.cancer.polyA.Basal.sample] %>% t,method='spearman') %>% c
Basal.cor.df <- rbind(
                  data.frame(cor.value = Basal.PRRX1.ECM.cor.vec,gene.type = 'ECM.gene'),
                  data.frame(cor.value = Basal.PRRX1.other.cor.vec,gene.type = 'other.gene')
                  )
Basal.cor.df$subtype <- 'Basal-like'

LumB.PRRX1.ECM.cor.vec  <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.LumB.sample] ,log2.fpkm.matrix[ECM.gene.ensemble.id,pure.TCGA.breast.cancer.polyA.LumB.sample] %>% t,method='spearman') %>% c
LumB.PRRX1.other.cor.vec <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.LumB.sample] ,log2.fpkm.matrix[other.gene,pure.TCGA.breast.cancer.polyA.LumB.sample] %>% t,method='spearman') %>% c
LumB.cor.df <- rbind(
  data.frame(cor.value = LumB.PRRX1.ECM.cor.vec,gene.type = 'ECM.gene'),
  data.frame(cor.value = LumB.PRRX1.other.cor.vec,gene.type = 'other.gene')
)

LumB.cor.df$subtype <- 'LuminalB'


Her2.PRRX1.ECM.cor.vec  <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.Her2.sample] ,log2.fpkm.matrix[ECM.gene.ensemble.id,pure.TCGA.breast.cancer.polyA.Her2.sample] %>% t,method='spearman') %>% c
Her2.PRRX1.other.cor.vec <- cor(log2.fpkm.matrix[PRRX1,pure.TCGA.breast.cancer.polyA.Her2.sample] ,log2.fpkm.matrix[other.gene,pure.TCGA.breast.cancer.polyA.Her2.sample] %>% t,method='spearman') %>% c

Her2.cor.df <- rbind(
  data.frame(cor.value = Her2.PRRX1.ECM.cor.vec,gene.type = 'ECM.gene'),
  data.frame(cor.value = Her2.PRRX1.other.cor.vec,gene.type = 'other.gene')
)
Her2.cor.df$subtype <- 'Her2-enriched'

save(file='client-side/output/PRRX1.R.output/PRRX1.RData',list=c('Basal.cor.df','LumB.cor.df','Her2.cor.df'))

draw.df <- rbind(Basal.cor.df,LumB.cor.df,Her2.cor.df)
draw.df$subtype <- factor(draw.df$subtype,levels = c('Basal-like','LuminalB','Her2-enriched'))
ggplot(draw.df,aes(x=subtype,y=cor.value,fill=gene.type)) + geom_violin(lwd=3)  + 
  scale_fill_manual(values = c('other.gene'='grey','ECM.gene'='blue')) + 
  theme_bw(base_size = 55) + theme(axis.title = element_text( size=25, face="bold"),
                                   axis.text.y  = element_text( size=55, face="bold"),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                   axis.line.x = element_line(colour = "black",size = 3),
                                   axis.line.y = element_line(colour = "black",size = 3),
                                   axis.text.x = element_text(angle = 45, hjust = 1,size=10, face="bold"),
                                   legend.position= 'none')            + xlab('')  + ylim (-1.2,1.2)

