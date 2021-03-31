source('client-side/code/DEBoost.R')

perform.DE.between.TWIST1.high.and.low <- function(){
  TWIST1        <- "ENSG00000026025" 
  TWIST1.expr   <- PRI.log2.tpm.matrix[TWIST1,] 
  q            <- quantile(TWIST1.expr)
  h.sample     <- names(TWIST1.expr)[TWIST1.expr >= q['75%']]
  l.sample     <- names(TWIST1.expr)[TWIST1.expr <= q['25%']]
  
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
TWIST1.DE.rs.list <- list()


load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.Basal.sample]
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.Basal.sample]
TWIST1.DE.rs.list[['BRCA.Basal']] <- perform.DE.between.TWIST1.high.and.low()


load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.Her2.sample]
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.Her2.sample]
TWIST1.DE.rs.list[['BRCA.Her2']] <- perform.DE.between.TWIST1.high.and.low()


load('server-side/RData//Breast Invasive Carcinoma.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output//Select.pure.sample.breast.cancer.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.breast.cancer.LumB.sample]
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.breast.cancer.LumB.sample]
TWIST1.DE.rs.list[['BRCA.Lumb']]  <- perform.DE.between.TWIST1.high.and.low()


load('client-side/output/Select.pure.sample.colorectal.cancer.R.output/Select.pure.sample.colorectal.cancer.RData')
load('server-side/RData/COLORECTAL_SRP029880.RData')
PRI.log2.read.count.matrix <- COLORECTAL_SRP029880_log2.read.count.matrix[,pure.PRI.colorectal.cancer.sample]
PRI.log2.tpm.matrix       <- COLORECTAL_SRP029880_log2.tpm.matrix[,pure.PRI.colorectal.cancer.sample]
TWIST1.DE.rs.list[['COAD']]       <- perform.DE.between.TWIST1.high.and.low()


load('client-side/output/Select.pure.sample.prostate.cancer.R.output/Select.pure.sample.prostate.cancer.RData')
load('server-side/RData//Prostate Adenocarcinoma.RData')
PRI.log2.read.count.matrix <- log2.read.count.matrix[,pure.PRI.prostate.cancer.sample]
PRI.log2.tpm.matrix       <- log2.tpm.matrix[,pure.PRI.prostate.cancer.sample]
TWIST1.DE.rs.list[['PRAD']]       <- perform.DE.between.TWIST1.high.and.low()


load('client-side/output/Select.pure.sample.NET.pancreatic.cancer.R.output/Select.pure.sample.NET.pancreatic.cancer.RData')
load('server-side/RData/GEP.NET.RData')
PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
PRI.log2.tpm.matrix       <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.pancreatic.cancer.sample]
TWIST1.DE.rs.list[['NET.PAAD']]   <- perform.DE.between.TWIST1.high.and.low()


load('client-side/output/Select.pure.sample.NET.si.cancer.R.output/Select.pure.sample.NET.si.cancer.RData')
load('server-side/RData/GEP.NET.RData')
PRI.log2.read.count.matrix <- GEP.NET.log2.read.count.matrix[,pure.PRI.NET.si.cancer.sample]
PRI.log2.tpm.matrix       <- GEP.NET.log2.tpm.matrix[,pure.PRI.NET.si.cancer.sample]
TWIST1.DE.rs.list[['NET.SI']]     <- perform.DE.between.TWIST1.high.and.low()


gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id


up.gene                  <- lapply(TWIST1.DE.rs.list,function(x) intersect(x$up.gene, mapping.df$ensembl_gene_id)    ) %>% unlist %>% unique
up.gene.matrix           <- matrix(0,nrow = length(up.gene),ncol = length(TWIST1.DE.rs.list))
rownames(up.gene.matrix) <- up.gene
colnames(up.gene.matrix) <- names(TWIST1.DE.rs.list)
for(i in 1:length(TWIST1.DE.rs.list) )  {
  gene   <- intersect(TWIST1.DE.rs.list[[i]]$up.gene, mapping.df$ensembl_gene_id)
  cancer <- names(TWIST1.DE.rs.list)[i]
  up.gene.matrix[gene,cancer] <- 1
}
up.gene.freq <- apply(up.gene.matrix,1,sum)
TWIST1.up.gene.freq <- up.gene.freq

dn.gene                  <- lapply(TWIST1.DE.rs.list,function(x) intersect(x$dn.gene, mapping.df$ensembl_gene_id)) %>% unlist %>% unique
dn.gene.matrix           <- matrix(0,nrow = length(dn.gene),ncol = length(TWIST1.DE.rs.list))
rownames(dn.gene.matrix) <- dn.gene
colnames(dn.gene.matrix) <- names(TWIST1.DE.rs.list)
for(i in 1:length(TWIST1.DE.rs.list) )  {
  gene   <- intersect(TWIST1.DE.rs.list[[i]]$dn.gene,mapping.df$ensembl_gene_id)
  cancer <- names(TWIST1.DE.rs.list)[i]
  dn.gene.matrix[gene,cancer] <- 1
}
dn.gene.freq <- apply(dn.gene.matrix,1,sum)
TWIST1.dn.gene.freq <- dn.gene.freq

TWIST1.dn.signature <- names(TWIST1.dn.gene.freq)[TWIST1.dn.gene.freq >= 4]
TWIST1.up.signature <- names(TWIST1.dn.gene.freq)[TWIST1.up.gene.freq >= 4]

