tmp <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt", stringsAsFactors=FALSE)
rownames(tmp) <- tmp$gene_id
tmp$gene_id  <- NULL
tmp$gene_name <- NULL
tmp$gene_type <- NULL
sc.expr.matrix <- as.matrix(tmp)


GSE75688.sample.meta <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/GSE75688_final_sample_information.txt", stringsAsFactors=FALSE)
flag <- GSE75688.sample.meta$index == 'Tumor' & GSE75688.sample.meta$index2 == 'Tumor' & GSE75688.sample.meta$index3 == 'Tumor' & GSE75688.sample.meta$type == 'SC'
sample.id <- GSE75688.sample.meta$sample[flag]
sc.expr.matrix <- sc.expr.matrix[,sample.id]
########### basal subtype ##################
get.patient <- function(x){
    strsplit(x=x,split = '_')[[1]] %>% head(1)       
}
primary.patient <- c('BC07','BC08','BC09','BC10','BC11')#
#primary.patient <- c('BC03')

patient.id <- sapply(colnames(sc.expr.matrix),get.patient)
p.matrix   <- sc.expr.matrix[,patient.id %in% primary.patient]


metastasis.patient <- c('BC07LN')
patient.id <- sapply(colnames(sc.expr.matrix),get.patient)
m.matrix   <- sc.expr.matrix[,patient.id %in% metastasis.patient]

rs <- foreach(g= rownames(sc.expr.matrix),.combine='rbind') %do% {
    p.value <- wilcox.test(p.matrix[g,],m.matrix[g,])$p.value 
    delta <- median(m.matrix[g,]) - median(p.matrix[g,])
    data.frame(effect.size=delta,p.value=p.value)
}

DOT1L <- 'ENSG00000104885'
HES4 <- 'ENSG00000188290'

p.value.vec <- sapply(1:nrow(p.matrix),function(i) wilcox.test(p.matrix[i,],m.matrix[i,])$p.value)
names(p.value.vec) <- rownames(sc.expr.matrix)
p.value.vec <- p.value.vec[is.na(p.value.vec) == FALSE]
p.value.vec <- sort(p.value.vec)

get.gene.id <- function(x){
  strsplit(x=x,split = '\\.')[[1]] %>% head(1)       
}

names(p.value.vec) <- sapply(names(p.value.vec), get.gene.id)
