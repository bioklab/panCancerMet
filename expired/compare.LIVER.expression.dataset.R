load('server-side//RData//Liver.RData')
GTex.liver.log2.fpkm.matrix       <- log2.fpkm.matrix
load('server-side/RData/SRP174668_Healthy.RData')
load('server-side/RData/SRP186450_Healthy.RData')
load('server-side/RData/SRP163252_Liver.RData')
load('server-side/RData/SRP068976_PairNor.RData')
load('server-side/RData/SRP174991_PairNor.RData')



protein.coding.gene.id <- read.table("~/Project/BreastCancerMetaPotenial/client-side/meta.data/protein.coding.gene.id.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
protein.coding.gene.id <- intersect(protein.coding.gene.id,rownames(GTex.liver.log2.fpkm.matrix))
g                      <- protein.coding.gene.id


SRP174668.liver.expr <- apply(SRP174668_log2.fpkm.matrix[g,],1,median)
SRP186450.liver.expr <- apply(SRP186450_log2.fpkm.matrix[g,],1,median) 
SRP163252.liver.expr <- apply(SRP163252_log2.fpkm.matrix[g,],1,median) 
SRP068976.liver.expr <- apply(SRP068976_log2.fpkm.matrix[g,],1,median) # 50 samples
SRP174991.liver.expr <- apply(SRP174991_log2.fpkm.matrix[g,],1,median) # 35 samples
GTex.liver.expr      <- apply(GTex.liver.log2.fpkm.matrix[g,],1,median) 


combined.liver.specific.gene.expr.matrix <- cbind(SRP174668.liver.expr,SRP186450.liver.expr,SRP163252.liver.expr,SRP068976.liver.expr,SRP174991.liver.expr,GTex.liver.expr)
pairs(combined.liver.specific.gene.expr.matrix)
cor.matrix       <- cor(combined.liver.specific.gene.expr.matrix,method='spearman')
diag(cor.matrix) <- 0
dataset.W        <- apply(cor.matrix,1,sum) %>% sort(decreasing = TRUE)

save(file = 'client-side/output/compare.LIVER.expression.dataset.R.output/compare.LIVER.expression.dataset.RData',list=c('combined.liver.specific.gene.expr.matrix','dataset.W'))
#################################
SAA1 <- 'ENSG00000173432'
DHODH <- 'ENSG00000102967'
