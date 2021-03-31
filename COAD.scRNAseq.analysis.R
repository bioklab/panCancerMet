require(dplyr)

load('client-side/Data/GSE97693/COLORECTAL_GSE97693_single.RData')
flag <- GSE97693_FPKM$X %in% c('1-Mar','2-Mar')
GSE97693_FPKM <- GSE97693_FPKM[!flag,]
rownames(GSE97693_FPKM) <- GSE97693_FPKM$X
GSE97693_FPKM$X <- NULL
GSE97693_FPKM.matrix <- apply(GSE97693_FPKM,2,as.numeric)
rownames(GSE97693_FPKM.matrix) <- rownames(GSE97693_FPKM)
flag   <- GSE97693_FPKM_METADATA$Sample_class %in% c(' Post-treatment Liver Metastasis','Liver Metastasis')
gsm.id <- GSE97693_FPKM_METADATA$GSM[flag] %>% as.character()


load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id

g.vec <- mapping.df[COAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.outlier.gene,'symbol']


GSE97693_FPKM.matrix['CD5L',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['SLC13A5',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['HABP2',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['F11',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['IGF2',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['MARCO',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['CPN1',gsm.id] %>% sort %>% boxplot
GSE97693_FPKM.matrix['INS-IGF2',gsm.id] %>% sort %>% boxplot


COAD.ectopic.liver.gene.expr.matrix <- GSE97693_FPKM.matrix[c('CD5L','SLC13A5','HABP2','F11','IGF2','MARCO','CPN1','INS-IGF2'),gsm.id]


GSE97693_FPKM.matrix['CPN1',gsm.id] %>% sort(decreasing = TRUE) %>% head
GSM2697027.expr.vec <- GSE97693_FPKM.matrix[,'GSM2697027'] # GSM2697027 has highest CPN1 expression



load('~/Project/Cancer2CellLine/server-side/RData/CCLE.RData')
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000   <- names(sort(rank.sd,decreasing =TRUE))[1:1000]

g            <- intersect(rownames(CCLE.log2.rpkm.matrix),mapping.df$ensembl_gene_id)
CC           <- CCLE.log2.rpkm.matrix[g,]
rownames(CC) <- mapping.df[g,'symbol']

m         <- mapping.df[CCLE.rna.seq.marker.gene.1000,'symbol']
m         <- m[is.na(m) == FALSE]
m         <- intersect(names(GSM2697027.expr.vec),m)
GSM2697027.cor.value <- cor(CC[m,],GSM2697027.expr.vec[m],method='spearman')

save(file='client-side/output/COAD.scRNAseq.analysis.R.output/COAD.scRNAseq.analysis.RData',list=c('GSM2697027.cor.value','COAD.ectopic.liver.gene.expr.matrix'))






# library(Matrix)
# matrix_dir = "/Users/keliu/Project/BreastCancerMetaPotenial/client-side/Data/GSE140312/filtered_feature_bc_matrix/"
# barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
# features.path <- paste0(matrix_dir, "features.tsv.gz")
# matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
# mat <- readMM(file = matrix.path)
# feature.names = read.delim(features.path, 
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
# barcode.names = read.delim(barcode.path, 
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
# colnames(mat) = barcode.names$V1
# rownames(mat) = feature.names$V1
# 
# 
# load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
# load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')
# 
# c.gene <- intersect(NET.PAAD.DE.rs$tumor.intrinsic.DE.gene.rs$up.outlier.gene, NET.SI.DE.rs$tumor.intrinsic.DE.gene.rs$up.outlier.gene)
# 
# r.matrix <- mat[c.gene,] %>% as.matrix



# library(dplyr)
# library(Seurat)
# library(patchwork)
# 
# # Load the PBMC dataset
# pbmc.data <- Read10X(data.dir = matrix_dir)
# # Initialize the Seurat object with the raw (non-normalized data).
# pbmc.counts <- Read10X(data.dir = matrix_dir)
# pbmc <- CreateSeuratObject(counts = pbmc.counts)
# pbmc <- NormalizeData(object = pbmc)
# pbmc <- FindVariableFeatures(object = pbmc)
# pbmc <- ScaleData(object = pbmc,features = rownames(pbmc))
# pbmc <- RunPCA(object = pbmc)
# pbmc <- FindNeighbors(object = pbmc)
# pbmc <- FindClusters(object = pbmc)
# pbmc <- RunTSNE(object = pbmc)
# DimPlot(object = pbmc, reduction = "tsne")
# 
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2




