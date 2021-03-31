
# Here we want to confirm maker genes are specificically expressed in the assocaited cell types
load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')
marker.gene <- c('ENSG00000167286', #CD3D, T cell
                 'ENSG00000010610', #CD4,  T cell
                 'ENSG00000153563', #CD8A, T cell
                 'ENSG00000172116', #CD8B, T cell   
                 'ENSG00000177455', #CD19, B cell
                 'ENSG00000156738', #CD20, B cell
                 'ENSG00000149294', #CD56, NK cell
                 'ENSG00000005961'  #CD41, Platlet
)
flag                  <- grepl(x=colnames(CCLE.log2.rpkm.matrix),pattern = 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE')
CCLE.log2.rpkm.matrix <- CCLE.log2.rpkm.matrix[,!flag]
apply(CCLE.log2.rpkm.matrix[marker.gene,],1,median)
boxplot(CCLE.log2.rpkm.matrix[marker.gene,] %>% t)

CCLE.log2.rpkm.matrix['ENSG00000177455',] %>% sort(decreasing = TRUE) %>% View
CCLE.log2.rpkm.matrix['ENSG00000005961',] %>% sort %>% View
