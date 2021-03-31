source('client-side/code/DEBoost.R')
load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
PRRX1        <- "ENSG00000116132" 

perform.DE.between.PRRX1.high.and.low <- function(){
    PRRX1        <- "ENSG00000116132" 
    PRRX1.expr   <- PRI.log2.tpm.matrix[PRRX1,] 
    q            <- quantile(PRRX1.expr)
    h.sample     <- names(PRRX1.expr)[PRRX1.expr >= q['75%']]
    l.sample     <- names(PRRX1.expr)[PRRX1.expr <= q['25%']]

    rs <- perform.DE.analysis.between.TRE.and.CON(CON.log2.tpm.matrix = PRI.log2.tpm.matrix[,h.sample],
                                                  TRE.log2.tpm.matrix = PRI.log2.tpm.matrix[,l.sample],
                                                  CON.log2.read.count.matrix = PRI.log2.read.count.matrix[,h.sample],
                                                  TRE.log2.read.count.matrix = PRI.log2.read.count.matrix[,l.sample]
    )

    up.gene <- rownames(rs)[rs$log2FoldChange >  1 & rs$padj < 0.05]
    dn.gene <- rownames(rs)[rs$log2FoldChange < -1 & rs$padj < 0.05]

    w.test.p.value <- foreach(g=c(up.gene,dn.gene),.combine='c') %do% {
        wilcox.test(PRI.log2.tpm.matrix[g,h.sample],PRI.log2.tpm.matrix[g,l.sample])$p.value  
    }
    names(w.test.p.value) <- c(up.gene,dn.gene)
    de.gene               <- names(w.test.p.value)[w.test.p.value < 0.05]
    up.gene               <- intersect(up.gene,de.gene)
    dn.gene               <- intersect(dn.gene,de.gene)
    list(up.gene=up.gene,dn.gene=dn.gene)
}
PRRX1.DE.rs.list <- list()
PRRX1.expr.rs.list <- list()


load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.Basal.sample]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]
#PRRX1.DE.rs.list[['BRCA.Basal']] <- perform.DE.between.PRRX1.high.and.low()

MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.Basal.sample]
df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
df$disease <- 'BRCA.Basal'
PRRX1.expr.rs.list[['BRCA.Basal']] <- df




load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.Her2.sample]
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.Her2.sample]
#PRRX1.DE.rs.list[['BRCA.Her2']] <- perform.DE.between.PRRX1.high.and.low()

MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.Her2.sample]
df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
df$disease <- 'BRCA.Her2'
PRRX1.expr.rs.list[['BRCA.Her2']] <- df

load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.LumB.sample]
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.LumB.sample]
#PRRX1.DE.rs.list[['BRCA.Lumb']]  <- perform.DE.between.PRRX1.high.and.low()

MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.MET.breast.cancer.LumB.sample]
df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
df$disease <- 'BRCA.LumB'
PRRX1.expr.rs.list[['BRCA.LumB']] <- df




load('client-side/output/Select.pure.sample.colorectal.cancer.R.output/Select.pure.sample.colorectal.cancer.RData')
load('server-side/RData/COLORECTAL_SRP029880.RData')
PRI.log2.read.count.matrix <- COLORECTAL_SRP029880_log2.read.count.matrix[,pure.PRI.colorectal.cancer.sample]
PRI.log2.tpm.matrix        <- COLORECTAL_SRP029880_log2.tpm.matrix[,pure.PRI.colorectal.cancer.sample]
#PRRX1.DE.rs.list[['COAD']]       <- perform.DE.between.PRRX1.high.and.low()


MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.MET.colorectal.cancer.sample]
df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
df$disease <- 'COAD'
PRRX1.expr.rs.list[['COAD']] <- df



load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')
load('server-side/RData//Prostate Adenocarcinoma.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.prostate.cancer.sample]
PRI.log2.tpm.matrix        <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]
#PRRX1.DE.rs.list[['PRAD']]       <- perform.DE.between.PRRX1.high.and.low()


MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.MET.prostate.cancer.sample]
df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
df$disease <- 'PRAD'
PRRX1.expr.rs.list[['PRAD']] <- df


load('client-side/output/Select.pure.sample.NET.pancreatic.cancer.R.output/Select.pure.sample.NET.pancreatic.cancer.RData')
load('server-side/RData/GEP.NET.RData')
PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
PRI.log2.tpm.matrix       <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
#PRRX1.DE.rs.list[['NET.PAAD']]   <- perform.DE.between.PRRX1.high.and.low()


MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
df$disease <- 'NET.PAAD'
PRRX1.expr.rs.list[['NET.PAAD']] <- df



load('client-side/output/Select.pure.sample.NET.si.cancer.R.output/Select.pure.sample.NET.si.cancer.RData')
load('server-side/RData/GEP.NET.RData')
PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.si.cancer.sample]
PRI.log2.tpm.matrix        <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.si.cancer.sample]
#PRRX1.DE.rs.list[['NET.SI']]     <- perform.DE.between.PRRX1.high.and.low()

MET.log2.tpm.matrix        <- MET500.log2.tpm.matrix[,pure.MET.NET.si.cancer.sample]
df         <- data.frame(expr = c(PRI.log2.tpm.matrix[PRRX1,],MET.log2.tpm.matrix[PRRX1,]), source = c(rep('PRI',PRI.log2.tpm.matrix %>% ncol),rep('MET',MET.log2.tpm.matrix %>% ncol)))  
df$disease <- 'NET.SI'
PRRX1.expr.rs.list[['NET.SI']] <- df



gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id


up.gene                  <- lapply(PRRX1.DE.rs.list,function(x) intersect(x$up.gene, mapping.df$ensembl_gene_id)    ) %>% unlist %>% unique
up.gene.matrix           <- matrix(0,nrow = length(up.gene),ncol = length(PRRX1.DE.rs.list))
rownames(up.gene.matrix) <- up.gene
colnames(up.gene.matrix) <- names(PRRX1.DE.rs.list)
for(i in 1:length(PRRX1.DE.rs.list) )  {
  gene   <- intersect(PRRX1.DE.rs.list[[i]]$up.gene, mapping.df$ensembl_gene_id)
  cancer <- names(PRRX1.DE.rs.list)[i]
  up.gene.matrix[gene,cancer] <- 1
}
up.gene.freq <- apply(up.gene.matrix,1,sum)
PRRX1.up.gene.freq <- up.gene.freq

dn.gene                  <- lapply(PRRX1.DE.rs.list,function(x) intersect(x$dn.gene, mapping.df$ensembl_gene_id)) %>% unlist %>% unique
dn.gene.matrix           <- matrix(0,nrow = length(dn.gene),ncol = length(PRRX1.DE.rs.list))
rownames(dn.gene.matrix) <- dn.gene
colnames(dn.gene.matrix) <- names(PRRX1.DE.rs.list)
for(i in 1:length(PRRX1.DE.rs.list) )  {
  gene   <- intersect(PRRX1.DE.rs.list[[i]]$dn.gene,mapping.df$ensembl_gene_id)
  cancer <- names(PRRX1.DE.rs.list)[i]
  dn.gene.matrix[gene,cancer] <- 1
}
dn.gene.freq <- apply(dn.gene.matrix,1,sum)
PRRX1.dn.gene.freq <- dn.gene.freq

PRRX1.dn.signature <- names(PRRX1.dn.gene.freq)[PRRX1.dn.gene.freq >= 6]
PRRX1.up.signature <- names(PRRX1.up.gene.freq)[PRRX1.up.gene.freq >= 6]

pooled.PRRX1.expr.df <- foreach(df = PRRX1.expr.rs.list,.combine='rbind') %do% {
  df
}

ggplot(pooled.PRRX1.expr.df,aes(x=disease,fill=source,y=expr)) + geom_boxplot()
