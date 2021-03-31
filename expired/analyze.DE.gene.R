#Single Cell paper:
## https://www.nature.com/articles/s41467-018-06052-0   Unravelling subclonal heterogeneity and aggressive disease states in TNBC through single-cell RNA-seq
## https://www.nature.com/articles/ncomms15081          Single-cell RNA-seq enables comprehensive tumour and immune cell profiling in primary breast cancer

##Landscape of tumor-infiltrating T cell repertoire of human cancers

# Let us check DE gene expr in primary tumor with high purity.
# If it is low, then we think DE is a reflection of TME

##Inhibition of histone methyltransferase DOT1L silences ERα gene and blocks proliferation of antiestrogen-resistant breast cancer cells
##In vivo screening identifies GATAD2B as a metastasis driver in KRAS-driven lung cancer. Nature communications
##Ref:FOXK1 and FOXK2 regulate aerobic glycolysis
##Ref:FOXK transcription factors: Regulation and critical role in cancer
##Ref:CDK11 Complexes Promote Pre-mRNA Splicing*
##Ref:Characterization of cyclin L1 and L2 interactions with CDK11 and splicing factors: influence of cyclin L isoforms on splice site selection.
##Ref:Cyclin-Dependent Kinase 11 (CDK11) Is Required for Ovarian Cancer Cell Growth In Vitro and In Vivo, and Its Inhibition Causes Apoptosis and Sensitizes Cells to Paclitaxel

##Ref:The Histone Methyltransferase DOT1L Promotes Neuroblastoma by Regulating Gene Transcription
##Ref:DOT1L cooperates with the c-Myc-p300 complex to epigenetically derepress CDH1 transcription factors in breast cancer progression
##Ref:Inhibition of histone methyltransferase DOT1L silences ERα gene and blocks proliferation of antiestrogen-resistant breast cancer cells
##


require(AnnotationDbi)
library(org.Hs.eg.db)
require(plyr)
require(dplyr)
require(stringr)
require(bedr)
require(CePa)
require(clusterProfiler)
###########################################################################################
### get a list of immune genes ###############
###########################################################################################

median.tpm.matrix <- read.gct(file='client-side/Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct')
get.gene.id <- function(x) {
    strsplit(x=x,split='\\.') %>% unlist %>% head(1)
}
rownames(median.tpm.matrix) <- sapply(rownames(median.tpm.matrix),get.gene.id)

is.immune.gene <- function(x) {
    y    <- log2(x+1)    
    tmp  <- quantile(y)
    upper <- tmp[4] + 1.5* IQR(y)
    outlier.tissue <- names(y)[ y > upper]  
    immune.tissue    <- c('Spleen','Whole.Blood')
    if(sum(outlier.tissue %in% immune.tissue) ==2){
        return(TRUE)
    }else{
        return(FALSE)  
    }
}
flag               <- apply(median.tpm.matrix,1,is.immune.gene)
immune.gene        <- rownames(median.tpm.matrix)[flag]
names(immune.gene) <- NULL
write.csv(x=immune.gene,file='client-side/output/analyze.DE.gene.R.output/immune.gene.list.csv',quote=FALSE)




###########################################################################################
#### remove immune genes, then check the overlappings of DE genes among different subtypes ###############
###########################################################################################

lumb.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv",  stringsAsFactors=FALSE)$x
basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
her2.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv",  stringsAsFactors=FALSE)$x
lumb.dn.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.dn.csv",  stringsAsFactors=FALSE)$x
her2.dn.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.dn.csv",  stringsAsFactors=FALSE)$x
basal.dn.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.dn.csv", stringsAsFactors=FALSE)$x

(intersect( c(lumb.up.gene,lumb.dn.gene), immune.gene)  %>% length ) / ( c(lumb.up.gene,lumb.dn.gene) %>% length )
(intersect( c(her2.up.gene,her2.dn.gene), immune.gene)  %>% length ) / ( c(her2.up.gene,her2.dn.gene) %>% length )
(intersect( c(basal.up.gene,basal.dn.gene), immune.gene)  %>% length ) / ( c(basal.up.gene,basal.dn.gene) %>% length )


basal.dn.gene    <- setdiff(basal.dn.gene,immune.gene)
basal.up.gene    <- setdiff(basal.up.gene,immune.gene)
her2.dn.gene     <- setdiff(her2.dn.gene, immune.gene)
her2.up.gene     <- setdiff(her2.up.gene, immune.gene)
lumb.dn.gene     <- setdiff(lumb.dn.gene, immune.gene)
lumb.up.gene     <- setdiff(lumb.up.gene, immune.gene)


write.csv(x=lumb.up.gene,  file = 'client-side/output/analyze.DE.gene.R.output/lumb.up.gene.immune.excluded.csv',  quote=FALSE)
write.csv(x=basal.up.gene, file = 'client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.csv', quote=FALSE)
write.csv(x=her2.up.gene,  file = 'client-side/output/analyze.DE.gene.R.output/her2.up.gene.immune.excluded.csv',  quote=FALSE)
write.csv(x=lumb.dn.gene,  file = 'client-side/output/analyze.DE.gene.R.output/lumb.dn.gene.immune.excluded.csv',  quote=FALSE)
write.csv(x=basal.dn.gene, file = 'client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.csv', quote=FALSE)
write.csv(x=her2.dn.gene,  file = 'client-side/output/analyze.DE.gene.R.output/her2.dn.gene.immune.excluded.csv',  quote=FALSE)



library(gplots)
venn(list(lumb=lumb.up.gene,her2=her2.up.gene,basal=basal.up.gene))
venn(list(lumb=lumb.dn.gene,her2=her2.dn.gene,basal=basal.dn.gene))


###########################################################################################
########  Convert ensemble ID to gene symbol, save the gene symbol list and then perform GO enrichment analysis  #############
###########################################################################################
lumb.up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = lumb.up.gene) 
basal.up.gene.annotation.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = basal.up.gene)
her2.up.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = her2.up.gene)
lumb.dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = lumb.dn.gene)
basal.dn.gene.annotation.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = basal.dn.gene)
her2.dn.gene.annotation.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = her2.dn.gene)

write.csv(x=lumb.up.gene.annotation.df, file = 'client-side/output/analyze.DE.gene.R.output/lumb.up.gene.immune.excluded.annotation.df.csv',  quote=FALSE)
write.csv(x=basal.up.gene.annotation.df,file = 'client-side/output/analyze.DE.gene.R.output/basal.up.gene.immune.excluded.annotation.df.csv', quote=FALSE)
write.csv(x=her2.up.gene.annotation.df,  file = 'client-side/output/analyze.DE.gene.R.output/her2.up.gene.immune.excluded.annotation.df.csv', quote=FALSE)

write.csv(x=lumb.dn.gene.annotation.df,  file = 'client-side/output/analyze.DE.gene.R.output/lumb.dn.gene.immune.excluded.annotation.df.csv', quote=FALSE)
write.csv(x=basal.dn.gene.annotation.df,file = 'client-side/output/analyze.DE.gene.R.output/basal.dn.gene.immune.excluded.annotation.df.csv', quote=FALSE)
write.csv(x=her2.dn.gene.annotation.df,  file = 'client-side/output/analyze.DE.gene.R.output/her2.dn.gene.immune.excluded.annotation.df.csv', quote=FALSE)



###########################################################################
#### GO enrichment 
###########################################################################
compute.jaccard.index  <- function(x,y){
  ( intersect(x,y) %>% length )/( unique(c(x,y))   %>% length )

}
get.gene.list <- function(x) {
    strsplit(x,split='\\/')[[1]] %>% unlist    
}



enriched.term.filtering <- function(x) {
    result    <- x %>% dplyr::filter(Count >=5 & qvalue <= 0.2 & pvalue<=0.001) %>% dplyr::arrange(pvalue) 
    if(nrow(result) == 0){
        return(NA)    
    }
    gene.str  <- result$geneID %>% unique
    idx.vec   <- sapply(gene.str,function(x) which(result$geneID == x) %>% min  )
    result    <- result[idx.vec,]
    gene.list <- lapply(result$geneID,get.gene.list)
    sim.matrix <- matrix(data=0,nrow=nrow(result),ncol=nrow(result))
    for(i in 1:nrow(sim.matrix)){
        for(j in 1:nrow(sim.matrix)){
            sim.matrix[i,j] <- compute.jaccard.index(gene.list[[i]],gene.list[[j]])    
        }
    }
    diag(sim.matrix) <- 0
    
    
    bit.vec   <- rep(1,nrow(result))
    for(i in 1:length(bit.vec)) {
        if(bit.vec[i] == 1){
            cover.index <- which(sim.matrix[i,] >= 0.4)
            bit.vec[cover.index] <- 0
        }        
    }
    result[bit.vec == 1,]
    
}

basal.up.gene.symbol  <- basal.up.gene.annotation.df$SYMBOL
basal.up.BP           <- enrichGO(gene=basal.up.gene.symbol,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
basal.up.BP.filtered  <- enriched.term.filtering(basal.up.BP)

basal.dn.gene.symbol  <- basal.dn.gene.annotation.df$SYMBOL
basal.dn.BP           <- enrichGO(gene=basal.dn.gene.symbol,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
basal.dn.BP.filtered  <- enriched.term.filtering(basal.dn.BP)


her2.up.gene.symbol  <- her2.up.gene.annotation.df$SYMBOL
her2.up.BP           <- enrichGO(gene=her2.up.gene.symbol,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
her2.up.BP.filtered  <- enriched.term.filtering(her2.up.BP)

her2.dh.gene.symbol  <- her2.dn.gene.annotation.df$SYMBOL
her2.dn.BP           <- enrichGO(gene=her2.dh.gene.symbol,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
her2.dn.BP.filtered  <- enriched.term.filtering(her2.dn.BP)



lumb.up.gene.symbol  <- lumb.up.gene.annotation.df$SYMBOL
lumb.up.BP           <- enrichGO(gene=lumb.up.gene.symbol,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
lumb.up.BP.filtered  <- enriched.term.filtering(lumb.up.BP)

lumb.dn.gene.symbol  <- lumb.dn.gene.annotation.df$SYMBOL
lumb.dn.BP           <- enrichGO(gene=lumb.dn.gene.symbol,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
lumb.dn.BP.filtered  <- enriched.term.filtering(lumb.dn.BP)




tmp                    <- c(lumb.dn.gene, her2.dn.gene,basal.dn.gene) %>% table %>% as.data.frame
c.gene                 <- tmp$.[tmp$Freq >= 2] %>% as.character()
c.gene.symbol          <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = c.gene)$SYMBOL
#c.gene.symbol          <- setdiff(c.gene.symbol,'MICOS10-NBL1') # well, ENSG00000158747 has two corresponding symbols
c.dn.gene.BP           <- enrichGO(gene=c.gene.symbol,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
c.dn.gene.BP.filtered  <- enriched.term.filtering(c.dn.gene.BP)
c.dn.gene              <- c.gene

tmp                    <- c(lumb.up.gene, her2.up.gene,basal.up.gene) %>% table %>% as.data.frame
c.gene                 <- tmp$.[tmp$Freq >= 2] %>% as.character()
c.gene.symbol          <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = c.gene)$SYMBOL
c.up.gene.BP           <- enrichGO(gene=c.gene.symbol,keyType='SYMBOL',OrgDb=org.Hs.eg.db,ont='BP',qvalueCutoff = 0.2,pvalueCutoff=0.001,maxGSSize = 200)@result
c.up.gene.BP.filtered  <- enriched.term.filtering(c.up.gene.BP)
c.up.gene              <- c.gene

save(file='client-side/output/analyze.DE.gene.R.output/common.DE.gene.RData',
     list=c('c.up.gene.BP.filtered','c.up.gene.BP','c.dn.gene.BP.filtered','c.dn.gene.BP','c.up.gene','c.dn.gene')
)



save(file='client-side/output/analyze.DE.gene.R.output/analyze.DE.gene.RData',list=c('basal.up.BP.filtered','basal.dn.BP.filtered','her2.up.BP.filtered','her2.dn.BP.filtered','lumb.up.BP.filtered','lumb.dn.BP.filtered',
                                                                                     'basal.up.BP',         'basal.dn.BP',         'her2.up.BP',         'her2.dn.BP',         'lumb.up.BP',         'lumb.dn.BP'
)
)



########################################################################################    
#####################  Trash code ###################################################### 
########################################################################################    











################################################################################################
## Check depmap score of DE genes
################################################################################################################

# load('client-side/output/organize.cancer.dependency.R.output/organize.cancer.dependency.RData')
# dep.score <- apply(Achilles_gene_effect.matrix,2,median)
# 
# basal.up.gene.dep.score <- dep.score[basal.up.gene.annotation.df$SYMBOL %>% as.character()]  %>% sort
# lumb.up.gene.dep.score  <- dep.score[lumb.up.gene.annotation.df$SYMBOL  %>% as.character()]  %>% sort
# her2.up.gene.dep.score  <- dep.score[her2.up.gene.annotation.df$SYMBOL  %>% as.character()]  %>% sort
# boxplot(basal.up.gene.dep.score,lumb.up.gene.dep.score,her2.up.gene.dep.score)
# 
# 
# basal.dn.gene.dep.score <- dep.score[basal.dn.gene.annotation.df$SYMBOL %>% as.character()]  %>% sort
# lumb.dn.gene.dep.score  <- dep.score[lumb.dn.gene.annotation.df$SYMBOL  %>% as.character()]  %>% sort
# her2.dn.gene.dep.score  <- dep.score[her2.dn.gene.annotation.df$SYMBOL  %>% as.character()]  %>% sort
# boxplot(basal.dn.gene.dep.score,lumb.dn.gene.dep.score,her2.dn.gene.dep.score)
# 
# boxplot(basal.up.gene.dep.score,basal.dn.gene.dep.score)
# boxplot(her2.up.gene.dep.score,her2.dn.gene.dep.score)
# boxplot(lumb.up.gene.dep.score,lumb.dn.gene.dep.score)
# 
# library(mixtools)
# x       <- dep.score[is.na(dep.score) == FALSE]
# mixmdl  <- normalmixEM(x)
# plot(mixmdl,which=2)

######## compute jaccard index ###########
# compute.jaccard.index  <- function(x,y){
#   ( intersect(x,y) %>% length )/( unique(c(x,y))   %>% length )
#   
# }
# 
# up.gene.jaccard.index <- c(compute.jaccard.index(basal.up.gene,her2.up.gene),
#                            compute.jaccard.index(basal.up.gene,lumb.up.gene),
#                            compute.jaccard.index(her2.up.gene,lumb.up.gene)
#                            )
# dn.gene.jaccard.index <- c(compute.jaccard.index(basal.dn.gene,her2.dn.gene),
#                            compute.jaccard.index(basal.dn.gene,lumb.dn.gene),
#                            compute.jaccard.index(her2.dn.gene,lumb.dn.gene)
# )



############## highly overlapped genes, DE in all four subtypes ###########

# tmp                 <- c(lumb.up.gene,basal.up.gene,her2.up.gene)
# freq.table          <- table(tmp) %>% as.data.frame
# high.up.ov.gene     <- freq.table$tmp[freq.table$Freq == 3] %>% as.character()
# high.up.ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = high.up.ov.gene)
# high.up.ov.gene.df  <- high.up.ov.gene.df[complete.cases(high.up.ov.gene.df),]
# 
# tmp                 <- c(lumb.dn.gene,basal.dn.gene,her2.dn.gene)
# freq.table          <- table(tmp) %>% as.data.frame
# high.dn.ov.gene     <- freq.table$tmp[freq.table$Freq == 3] %>% as.character()
# high.dn.ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = high.dn.ov.gene)
# high.dn.ov.gene.df  <- high.dn.ov.gene.df[complete.cases(high.dn.ov.gene.df),]


# ##############  overlapped genes, DE in all four subtypes ###########
# 
# tmp            <- c(luma.up.gene,lumb.up.gene,basal.up.gene,her2.up.gene)
# freq.table     <- table(tmp) %>% as.data.frame
# up.ov.gene     <- freq.table$tmp[freq.table$Freq >= 2] %>% as.character()
# up.ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.ov.gene)
# up.ov.gene.df  <- up.ov.gene.df[complete.cases(up.ov.gene.df),]
# 
# tmp            <- c(luma.dn.gene,lumb.dn.gene,basal.dn.gene,her2.dn.gene)
# freq.table     <- table(tmp) %>% as.data.frame
# dn.ov.gene     <- freq.table$tmp[freq.table$Freq >= 2] %>% as.character()
# dn.ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = dn.ov.gene)
# dn.ov.gene.df  <- dn.ov.gene.df[complete.cases(dn.ov.gene.df),]




###################################################################################################
##### Clustering of subtypes based on DE genes
###################################################################################################
# tmp            <- c(lumb.up.gene,basal.up.gene,her2.up.gene)
# freq.table     <- table(tmp) %>% as.data.frame
# up.ov.gene     <- freq.table$tmp[freq.table$Freq >= 2] %>% as.character()
# up.ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = up.ov.gene)
# up.ov.gene.df  <- up.ov.gene.df[complete.cases(up.ov.gene.df),]
# c.TF           <- intersect(up.ov.gene.df$SYMBOL,TF.list)
# 
# up.ov.gene.matrix <- foreach(g=up.ov.gene,.combine='rbind') %do% {
#     c(ifelse(g %in% luma.up.gene,1,0), 
#       ifelse(g %in% lumb.up.gene,1,0), 
#       ifelse(g %in% her2.up.gene,1,0), 
#       ifelse(g %in% basal.up.gene,1,0)
#     )
# }
# rownames(up.ov.gene.matrix) <- up.ov.gene
# colnames(up.ov.gene.matrix) <- c('LuminalA','LuminalB','Her2-enriched','Basal')
# 
# tmp          <- c(luma.dn.gene,lumb.dn.gene,basal.dn.gene,her2.dn.gene)
# freq.table   <- table(tmp) %>% as.data.frame
# dn.ov.gene   <- freq.table$tmp[freq.table$Freq >= 2] %>% as.character()
# dn.ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = dn.ov.gene)
# dn.ov.gene.df  <- dn.ov.gene.df[complete.cases(dn.ov.gene.df),]
# 
# dn.ov.gene.matrix <- foreach(g=dn.ov.gene,.combine='rbind') %do% {
#     c(ifelse(g %in% luma.dn.gene,1,0), 
#       ifelse(g %in% lumb.dn.gene,1,0), 
#       ifelse(g %in% her2.dn.gene,1,0), 
#       ifelse(g %in% basal.dn.gene,1,0)
#     )
# }
# rownames(dn.ov.gene.matrix) <- dn.ov.gene
# colnames(dn.ov.gene.matrix) <- c('LuminalA','LuminalB','Her2-enriched','Basal')
# 
# require(pheatmap)
# pheatmap(rbind(up.ov.gene.matrix,dn.ov.gene.matrix))
# 
# 
# luma.gene <- c(luma.up.gene,luma.dn.gene)
# lumb.gene <- c(lumb.up.gene,lumb.dn.gene)
# her2.gene <- c(her2.up.gene,her2.dn.gene)
# basal.gene <- c(basal.up.gene,basal.dn.gene)
# 
# compute.jaccard.index  <- function(x,y){
#    ( intersect(x,y) %>% length )/( unique(c(x,y))   %>% length )
#   
# }
# 
# jaccard.matrix <- matrix(data = 1,nrow=4,ncol=4)
# rownames(jaccard.matrix) <- c('LuminalA','LuminalB','Her2-enriched','Basal-like')
# colnames(jaccard.matrix) <- c('LuminalA','LuminalB','Her2-enriched','Basal-like')
# 
# 
# 
# jaccard.matrix['LuminalA','LuminalB']      <- compute.jaccard.index(luma.up.gene,lumb.up.gene)
# jaccard.matrix['LuminalA','Her2-enriched'] <- compute.jaccard.index(luma.up.gene,her2.up.gene)
# jaccard.matrix['LuminalA','Basal-like']         <- compute.jaccard.index(luma.up.gene,basal.up.gene)
# jaccard.matrix['LuminalB','Her2-enriched'] <- compute.jaccard.index(lumb.up.gene,her2.up.gene)
# jaccard.matrix['LuminalB','Basal-like']         <- compute.jaccard.index(lumb.up.gene,basal.up.gene)
# jaccard.matrix['Her2-enriched','Basal-like']    <- compute.jaccard.index(her2.up.gene,basal.up.gene)
# 
# 
# 
# 
# jaccard.matrix['LuminalB','LuminalA']      <-compute.jaccard.index(luma.dn.gene,lumb.dn.gene)
# jaccard.matrix['Her2-enriched','LuminalA'] <-compute.jaccard.index(luma.dn.gene,her2.dn.gene)
# jaccard.matrix['Basal-like','LuminalA']         <-compute.jaccard.index(luma.dn.gene,basal.dn.gene)
# jaccard.matrix['Her2-enriched','LuminalB'] <-compute.jaccard.index(lumb.dn.gene,her2.dn.gene)
# jaccard.matrix['Basal-like','LuminalB']         <-compute.jaccard.index(lumb.dn.gene,basal.dn.gene)
# jaccard.matrix['Basal-like','Her2-enriched']    <-compute.jaccard.index(her2.dn.gene,basal.dn.gene)
# 
# 
# up.gene.jaccard.matrix                   <- jaccard.matrix
# up.gene.jaccard.matrix[lower.tri(up.gene.jaccard.matrix)]        <- 0
# up.gene.jaccard.matrix                   <- up.gene.jaccard.matrix + t(up.gene.jaccard.matrix)
# diag(up.gene.jaccard.matrix)             <- 1
# rownames(up.gene.jaccard.matrix) <- c('LuminalA','LuminalB','Her2-enriched','Basal-like')
# colnames(up.gene.jaccard.matrix) <- c('LuminalA','LuminalB','Her2-enriched','Basal-like')
# up.gene.hclust.rs <- hclust(as.dist(1-up.gene.jaccard.matrix),method='average')
# plot(up.gene.hclust.rs,font=2,cex=2,cex.x=2)
# 
# 
# dn.gene.jaccard.matrix                   <- jaccard.matrix
# dn.gene.jaccard.matrix[upper.tri(dn.gene.jaccard.matrix)]        <- 0
# dn.gene.jaccard.matrix                   <- dn.gene.jaccard.matrix + t(dn.gene.jaccard.matrix)
# diag(dn.gene.jaccard.matrix)             <- 1
# rownames(dn.gene.jaccard.matrix) <- c('LuminalA','LuminalB','Her2-enriched','Basal-like')
# colnames(dn.gene.jaccard.matrix) <- c('LuminalA','LuminalB','Her2-enriched','Basal-like')
# dn.gene.hclust.rs <- hclust(as.dist(1-dn.gene.jaccard.matrix),,method='average')
# plot(dn.gene.hclust.rs,font=2,cex=1.5,xlab='')
# 
# 
# 
# ############## subtype specific up and dn regulated transcription regulators ########
# 
# basal.specific.up.TF <- setdiff(intersect(basal.up.gene.annotation.df$SYMBOL,TF.list),up.ov.gene.df$SYMBOL)
# lumb.specific.up.TF  <- setdiff(intersect(lumb.up.gene.annotation.df$SYMBOL,TF.list),up.ov.gene.df$SYMBOL)
# luma.specific.up.TF  <- setdiff(intersect(luma.up.gene.annotation.df$SYMBOL,TF.list),up.ov.gene.df$SYMBOL)
# her2.specific.up.TF  <- setdiff(intersect(her2.up.gene.annotation.df$SYMBOL,TF.list),up.ov.gene.df$SYMBOL)
# 
# basal.specific.dn.TF <- setdiff(intersect(basal.dn.gene.annotation.df$SYMBOL,TF.list),dn.ov.gene.df$SYMBOL)
# lumb.specific.dn.TF  <- setdiff(intersect(lumb.dn.gene.annotation.df$SYMBOL,TF.list),dn.ov.gene.df$SYMBOL)
# luma.specific.dn.TF  <- setdiff(intersect(luma.dn.gene.annotation.df$SYMBOL,TF.list),dn.ov.gene.df$SYMBOL)
# her2.specific.dn.TF  <- setdiff(intersect(her2.dn.gene.annotation.df$SYMBOL,TF.list),dn.ov.gene.df$SYMBOL)






#########
# targetable.gene <- read.table("~/Project/BreastCancerMetaPotenial/client-side/Data/targetable_GeneNames_chembl25.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)$V1
# basal.up.gene.dep.score[intersect(names(basal.up.gene.dep.score),targetable.gene)]
# lumb.up.gene.dep.score[intersect(names(lumb.up.gene.dep.score),targetable.gene)] %>% sort() # lumb, emphsize CDK11A. https://stm.sciencemag.org/content/11/509/eaaw8412
# her2.up.gene.dep.score[intersect(names(her2.up.gene.dep.score),targetable.gene)] %>% sort() # 

### should claim that higher expression, higher depedency?########



#######
# require(GSA)
# msigdb          <-  GSA.read.gmt("~/Project/Cancer2CellLine/client-side/meta.data/h.all.v6.1.entrez.gmt")
# genesets        <-  msigdb$genesets
# names(genesets) <-  msigdb$geneset.names 
# 
# EMT.gene         <- genesets$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
# EMT.gene.symbol  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENTREZID',columns = c('ENTREZID','SYMBOL'),keys = EMT.gene)$SYMBOL
# x                <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'GO',columns = c('GO','SYMBOL'),keys = c('GO:0010718'))$SYMBOL #GO:0010718:  positive regulation of epithelial to mesenchymal transition
# EMT.gene.symbol  <- c(EMT.gene.symbol,x) %>% unique
# x                <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'GO',columns = c('GO','SYMBOL'),keys = c('GO:0010719'))$SYMBOL #GO:0010719:  negative regulation of epithelial to mesenchymal transition
# EMT.gene.symbol  <- c(EMT.gene.symbol,x) %>% unique
# 
# 
# intersect(basal.up.gene.annotation.df$SYMBOL,EMT.gene.symbol)
# intersect(basal.dn.gene.annotation.df$SYMBOL,EMT.gene.symbol)
# 
# intersect(luma.up.gene.annotation.df$SYMBOL,EMT.gene.symbol)
# intersect(luma.dn.gene.annotation.df$SYMBOL,EMT.gene.symbol)
# 
# intersect(lumb.up.gene.annotation.df$SYMBOL,EMT.gene.symbol)
# intersect(lumb.dn.gene.annotation.df$SYMBOL,EMT.gene.symbol)
# 
# intersect(her2.up.gene.annotation.df$SYMBOL,EMT.gene.symbol)
# intersect(her2.dn.gene.annotation.df$SYMBOL,EMT.gene.symbol)
# 
# 
# 
# #ACTA2 mesemchymal marker genes
# 
# SNAI1 <- 'ENSG00000124216'
# SNAI2 <- 'ENSG00000019549'
# ZEB1  <- 'ENSG00000148516'
# ZEB2  <- 'ENSG00000169554'
# TWIST <- 'ENSG00000122691'
# CDH1  <- 'ENSG00000039068'
# EMT.TF <- c(SNAI1,SNAI2,ZEB1,ZEB2,TWIST,CDH1)
# 
# load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
# de.res.metastasis.liver.vs.breast.basal[EMT.TF,]
# de.res.liver.vs.breast.basal[EMT.TF,]
# 
# de.res.metastasis.liver.vs.breast.her2[EMT.TF,]
# de.res.liver.vs.breast.her2[EMT.TF,]
# 
# de.res.metastasis.liver.vs.breast.lumb[EMT.TF,]
# de.res.liver.vs.breast.lumb[EMT.TF,]
# 
# de.res.metastasis.liver.vs.breast.luma[EMT.TF,]
# de.res.liver.vs.breast.luma[EMT.TF,]
# 













#############################################################################################
## Assemble TF list
###########################################################################################
# file.name <- "client-side/Data/JASPAR.txt"
# conn      <- file(file.name,open="r")
# line      <- readLines(conn)
# get.TF.name <- function(x) {
#     if(grepl(x = x,pattern='>')){
#         l <- strsplit(x=x,split='\t')  %>% unlist
#         toupper(l[2])
#     }else{
#         NA  
#     }  
# }
# 
# TF.list <- sapply(line,get.TF.name)
# TF.list <- TF.list[is.na(TF.list) == FALSE]
# names(TF.list) <- NULL
# JASPAR.TF.list <- TF.list
# close(conn)
# 
# 
# file.name <- "client-side/Data/CISTROME.factor.line.txt"
# conn      <- file(file.name,open="r")
# line      <- readLines(conn)
# tmp       <- str_extract_all(pattern="id=\"[:alnum:]+\"",string=line[1],simplify = TRUE)
# get.TF.name <- function(x) {
#     x <- str_remove_all(string = x,pattern="\"")
#     x <- str_remove_all(string = x,pattern="id=")
#     x
# }
# TF.list <- sapply(tmp,get.TF.name)
# TF.list <- TF.list[is.na(TF.list) == FALSE]
# names(TF.list) <- NULL
# CISTROME.TF.list <- TF.list
# close(conn)
# 
# TF.list <- c(CISTROME.TF.list,JASPAR.TF.list) %>% unique



# ########## for dn-regulated genes
# luma.dn.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = luma.dn.gene)
# lumb.dn.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = lumb.dn.gene)
# basal.dn.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = basal.dn.gene)
# her2.dn.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = her2.dn.gene)
# 
# tmp         <- c(luma.dn.gene,lumb.dn.gene,basal.dn.gene,her2.dn.gene)
# freq.table  <- table(tmp) %>% as.data.frame
# ov.gene     <- freq.table$tmp[freq.table$Freq >= 2] %>% as.character()
# ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = ov.gene)
# ov.gene.df  <- ov.gene.df[complete.cases(ov.gene.df),]
# c.TF        <- intersect(ov.gene.df$SYMBOL,TF.list)
# 
# high.ov.gene     <- freq.table$tmp[freq.table$Freq >= 4] %>% as.character()
# high.ov.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = high.ov.gene)
# high.ov.gene.df  <- high.ov.gene.df[complete.cases(high.ov.gene.df),]
# 
# 
# basal.specific.TF <- setdiff(intersect(basal.dn.gene.df$SYMBOL,TF.list),ov.gene.df$SYMBOL)
# lumb.specific.TF  <- setdiff(intersect(lumb.dn.gene.df$SYMBOL,TF.list),ov.gene.df$SYMBOL)
# luma.specific.TF  <- setdiff(intersect(luma.dn.gene.df$SYMBOL,TF.list),ov.gene.df$SYMBOL)
# her2.specific.TF  <- setdiff(intersect(her2.dn.gene.df$SYMBOL,TF.list),ov.gene.df$SYMBOL)
# 
# 
# ####################
# load('~/Project/GeneFishing/RData/normalized.GTex.with.ncRNA.median.larger.than.0.1.RData')
# 
# breast.data   <- normalized.GTex.data[['Breast - Mammary Tissue']]
# g.vec         <- intersect(colnames(breast.data),basal.up.gene)
# outlier.vec   <- foreach(g = g.vec,.combine='c') %do% {
#     v              <- cor(breast.data[,g],breast.data,method='spearman') %>% c  
#     up.outlier     <- quantile(v)[4] + 1.5 * IQR(v)
#     up.outlier
# }
# names(outlier.vec) <- g.vec
# 
# cor.matrix <- cor(breast.data[,g.vec],method='spearman')
# 
# cm <- combn(g.vec,m = 2)
# 
# f.vec <- foreach(j=1:ncol(cm),.combine='c') %do% {
#     g1    <- cm[1,j]
#     g2    <- cm[2,j]
#     value <- cor.matrix[g1,g2]
#     ifelse(value > outlier.vec[g1] & value > outlier.vec[g2],TRUE,FALSE)
# }
# edge <- cm[,f.vec]
# 
# 
# edge[1,] <- basal.up.gene.df[edge[1,],'SYMBOL']
# edge[2,] <- basal.up.gene.df[edge[2,],'SYMBOL']
# 
# flag <- edge[1,] %in% TF.list | edge[2,] %in% TF.list
# t.gene <- c(edge[,flag]) %>% unique
# edge[,flag] %>% View
# 
# ############### load the biogrid network ################
# require(data.table)
# data <- fread(input = 'client-side/Data/BIOGRID-ALL-3.5.175.tab2.txt')
# data <- as.data.frame(data)
# human.network <- data[data$`Organism Interactor A` == 9606,]
# human.network <- human.network[human.network$`Entrez Gene Interactor A` != human.network$`Entrez Gene Interactor B`,] # remove self-self interactions

# 
# 
# load('client-side/output//TCGA.breast.cancer.meta.R.output//TCGA.breast.cancer.meta.RData')
# load('server-side/RData//Breast Invasive Carcinoma.RData')
# load('client-side/output/tumor.purity.based.on.cell.line.R.output/tumor.purity.based.on.cell.line.RData')
# 
# ####################### Basal subtype #######################
# df      <- read.csv('client-side/output/DE.breast.cancer.R.output//basal.up.csv',header = TRUE)
# up.gene <- df$x %>% as.character
# df      <- read.csv('client-side/output/DE.breast.cancer.R.output//basal.dn.csv',header = TRUE)
# dn.gene <- df$x %>% as.character
# up.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=up.gene, columns=c('ENTREZID','SYMBOL'), keytype="ENSEMBL")
# up.gene.annotation <- up.gene.annotation[complete.cases(up.gene.annotation),]
# dn.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=dn.gene, columns=c('ENTREZID','SYMBOL'), keytype="ENSEMBL")
# dn.gene.annotation <- dn.gene.annotation[complete.cases(dn.gene.annotation),]
# 
# gene <- dn.gene.annotation$ENTREZID %>% unique
# gene <- gene[is.na(gene) == FALSE]
# 
# flag <- (human.network$`Entrez Gene Interactor A` %in% gene) & (human.network$`Entrez Gene Interactor B` %in% gene)
# sub.network <- human.network[flag,8:9] %>% unique
# write.csv(x=sub.network,file = 'client-side/output/analyze.DE.gene.R.output/Basal.network.csv',quote = FALSE,row.names=FALSE)

# 
# ####################### LumB subtype #######################
# df      <- read.csv('client-side/output/DE.breast.cancer.R.output//lumb.up.csv',header = TRUE)
# up.gene <- df$x %>% as.character
# df      <- read.csv('client-side/output/DE.breast.cancer.R.output//lumb.dn.csv',header = TRUE)
# dn.gene <- df$x %>% as.character
# up.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=up.gene, columns=c('ENTREZID','SYMBOL'), keytype="ENSEMBL")
# dn.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=dn.gene, columns=c('ENTREZID','SYMBOL'), keytype="ENSEMBL")
# 
# gene <- up.gene.annotation$ENTREZID %>% unique
# gene <- gene[is.na(gene) == FALSE]
# 
# flag <- (human.network$`Entrez Gene Interactor A` %in% gene) & (human.network$`Entrez Gene Interactor B` %in% gene)
# sub.network <- human.network[flag,8:9] %>% unique
# write.csv(x=sub.network,file = 'client-side/output/analyze.DE.gene.R.output/LumB.network.csv',quote = FALSE,row.names=FALSE)

# 
# ####################### LumA subtype #######################
# df      <- read.csv('client-side/output/DE.breast.cancer.R.output//luma.up.csv',header = TRUE)
# up.gene <- df$x %>% as.character
# df      <- read.csv('client-side/output/DE.breast.cancer.R.output//luma.dn.csv',header = TRUE)
# dn.gene <- df$x %>% as.character
# up.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=up.gene, columns=c('ENTREZID','SYMBOL'), keytype="ENSEMBL")
# 
# gene <- up.gene.annotation$ENTREZID %>% unique
# gene <- gene[is.na(gene) == FALSE]
# 
# flag <- (human.network$`Entrez Gene Interactor A` %in% gene) & (human.network$`Entrez Gene Interactor B` %in% gene)
# sub.network <- human.network[flag,8:9] %>% unique
# write.csv(x=sub.network,file = 'client-side/output/analyze.DE.gene.R.output/LumA.network.csv',quote = FALSE,row.names=FALSE)

# 
# ####################### Her2 subtype #######################
# df      <- read.csv('client-side/output/DE.breast.cancer.R.output//her2.up.csv',header = TRUE)
# up.gene <- df$x %>% as.character
# df      <- read.csv('client-side/output/DE.breast.cancer.R.output//her2.dn.csv',header = TRUE)
# dn.gene <- df$x %>% as.character
# up.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=up.gene, columns=c('ENTREZID','SYMBOL'), keytype="ENSEMBL")
# 
# gene <- up.gene.annotation$ENTREZID %>% unique
# gene <- gene[is.na(gene) == FALSE]
# 
# flag <- (human.network$`Entrez Gene Interactor A` %in% gene) & (human.network$`Entrez Gene Interactor B` %in% gene)
# sub.network <- human.network[flag,8:9] %>% unique
# write.csv(x=sub.network,file = 'client-side/output/analyze.DE.gene.R.output/Her2.network.csv',quote = FALSE,row.names=FALSE)

# ###################
# 
# luma.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/luma.up.csv",  stringsAsFactors=FALSE)$x
# lumb.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/lumb.up.csv",  stringsAsFactors=FALSE)$x
# basal.up.gene    <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/basal.up.csv", stringsAsFactors=FALSE)$x
# her2.up.gene     <- read.csv("~/Project/BreastCancerMetaPotenial/client-side/output/DE.breast.cancer.R.output/her2.up.csv",  stringsAsFactors=FALSE)$x
# 
# luma.up.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = luma.up.gene)
# lumb.up.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = lumb.up.gene)
# basal.up.gene.df <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = basal.up.gene)
# her2.up.gene.df  <- AnnotationDbi::select(x = org.Hs.eg.db,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'),keys = her2.up.gene)
# 
# 
# write.csv(x=luma.up.gene.df,file = 'client-side/output/analyze.DE.gene.R.output/luma.up.gene.symbol.df.csv',quote=FALSE)
# write.csv(x=lumb.up.gene.df,file = 'client-side/output/analyze.DE.gene.R.output/lumb.up.gene.symbol.df.csv',quote=FALSE)
# write.csv(x=basal.up.gene.df,file = 'client-side/output/analyze.DE.gene.R.output/basal.up.gene.symbol.df.csv',quote=FALSE)
# write.csv(x=her2.up.gene.df,file = 'client-side/output/analyze.DE.gene.R.output/her2.up.gene.symbol.df.csv',quote=FALSE)
# 
# 
# 
# #dn.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=dn.gene, columns=c('ENTREZID','SYMBOL'), keytype="ENSEMBL")
# 
# # basal.up.gene <- up.gene
# # basal.dn.gene <- dn.gene
# # # Gene, encoded protein nucleus-located
# # # GO.nucleus <- 'GO:0005634'
# # # nucleus.gene.symbol <- up.gene.annotation$SYMBOL[up.gene.annotation$GOALL == GO.nucleus] %>% unique
# # # nucleus.gene.symbol <- nucleus.gene[is.na(nucleus.gene) == FALSE]
# 
# 
# 
# ######################## LumB subtype #######################
# df      <- read.csv('client-side/output/DE.R.output/lumb.up.csv',header = TRUE)
# up.gene <- df$x %>% as.character
# df      <- read.csv('client-side/output/DE.R.output/lumb.dn.csv',header = TRUE)
# dn.gene <- df$x %>% as.character
# 
# up.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=up.gene, columns=c('ENTREZID','SYMBOL','GOALL'), keytype="ENSEMBL")
# dn.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=dn.gene, columns=c('ENTREZID','SYMBOL','GOALL'), keytype="ENSEMBL")
# 
# 
# # Gene, encoded protein nucleus-located
# # GO.nucleus <- 'GO:0005634'
# # nucleus.gene.symbol <- up.gene.annotation$SYMBOL[up.gene.annotation$GOALL == GO.nucleus] %>% unique
# # nucleus.gene.symbol <- nucleus.gene.symbol[is.na(nucleus.gene.symbol) == FALSE]
# 
# lumb.up.gene <- up.gene
# lumb.dn.gene <- dn.gene
# 
# gene <- up.gene.annotation$ENTREZID %>% unique
# gene <- gene[is.na(gene) == FALSE]
# flag <- (human.network$`Entrez Gene Interactor A` %in% gene) & (human.network$`Entrez Gene Interactor B` %in% gene)
# sub.network <- human.network[flag,8:9] %>% unique
# 
# write.csv(x=sub.network,file = 'client-side/output/central.dogma.R.output/LumB.network.csv',quote = FALSE,row.names=FALSE)
# 
# ###################### Her2 subtype ######################
# df      <- read.csv('client-side/output/DE.R.output/her2.up.csv',header = TRUE)
# up.gene <- df$x %>% as.character
# df      <- read.csv('client-side/output/DE.R.output/her2.dn.csv',header = TRUE)
# dn.gene <- df$x %>% as.character
# 
# up.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=up.gene, columns=c('ENTREZID','SYMBOL','GOALL'), keytype="ENSEMBL")
# dn.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=dn.gene, columns=c('ENTREZID','SYMBOL','GOALL'), keytype="ENSEMBL")
# 
# 
# # Gene, encoded protein nucleus-located
# GO.nucleus <- 'GO:0005634'
# nucleus.gene.symbol <- up.gene.annotation$SYMBOL[up.gene.annotation$GOALL == GO.nucleus] %>% unique
# nucleus.gene.symbol <- nucleus.gene[is.na(nucleus.gene) == FALSE]
# 
# 
# 
# nucleus.gene <- up.gene.annotation$ENTREZID[up.gene.annotation$GOALL == GO.nucleus] %>% unique
# nucleus.gene <- nucleus.gene[is.na(nucleus.gene) == FALSE]
# 
# flag <- (human.network$`Entrez Gene Interactor A` %in% nucleus.gene) & (human.network$`Entrez Gene Interactor B` %in% nucleus.gene)
# sub.network <- human.network[flag,8:9] %>% unique
# 
# write.csv(x=sub.network,file = 'client-side/output/central.dogma.R.output/Her2.network.csv',quote = FALSE,row.names=FALSE)
# 
# her2.up.gene <- up.gene
# her2.dn.gene <- dn.gene
# 
# ###################### Luma subtype ############################
# df      <- read.csv('client-side/output/DE.R.output/luma.up.csv',header = TRUE)
# up.gene <- df$x %>% as.character
# df      <- read.csv('client-side/output/DE.R.output/luma.dn.csv',header = TRUE)
# dn.gene <- df$x %>% as.character
# 
# up.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=up.gene, columns=c('ENTREZID','SYMBOL','GOALL'), keytype="ENSEMBL")
# dn.gene.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys=dn.gene, columns=c('ENTREZID','SYMBOL','GOALL'), keytype="ENSEMBL")
# 
# 
# # Gene, encoded protein nucleus-located
# GO.nucleus <- 'GO:0005634'
# nucleus.gene.symbol <- up.gene.annotation$SYMBOL[up.gene.annotation$GOALL == GO.nucleus] %>% unique
# nucleus.gene.symbol <- nucleus.gene[is.na(nucleus.gene) == FALSE]
# 
# 
# 
# nucleus.gene <- up.gene.annotation$ENTREZID[up.gene.annotation$GOALL == GO.nucleus] %>% unique
# nucleus.gene <- nucleus.gene[is.na(nucleus.gene) == FALSE]
# 
# flag <- (human.network$`Entrez Gene Interactor A` %in% nucleus.gene) & (human.network$`Entrez Gene Interactor B` %in% nucleus.gene)
# sub.network <- human.network[flag,8:9] %>% unique
# 
# write.csv(x=sub.network,file = 'client-side/output/central.dogma.R.output/LumA.network.csv',quote = FALSE,row.names=FALSE)
# 
# luma.up.gene <- up.gene
# luma.dn.gene <- dn.gene


# subtype.sample <- TCGA.breast.cancer.polyA.Basal.sample
# expr.data      <- log2.fpkm.matrix[,subtype.sample]
# purity.vec     <- tumor.purity.based.on.cell.line.vec[subtype.sample] 
# 
# g <- c(up.gene,dn.gene)
# plot(x=purity.vec,y=expr.data['ENSG00000153563',])
# 
# cor.vec <- cor(expr.data[g,] %>% t,purity.vec,method='spearman')
# cor.vec <- c(cor.vec)
# names(cor.vec) <- g
# #cor.vec <- sort(cor.vec)
# 
# plot(x=purity.vec,y=expr.data['ENSG00000124875',])
# 
# boxplot(cor.vec[up.gene],cor.vec[dn.gene])
# 
# e.g <- names(cor.vec)[cor.vec > -0.3]
# 
# u.e.g <- names(cor.vec)[cor.vec < -0.3]
# 
# lumb.u.e.g <- u.e.g
# basal.u.e.g <- u.e.g
# 
# 
# 
# 
# 
# 
# tmp            <- as.integer(0.2 * length(purity.vec))
# purity.vec    <- purity.vec[1:tmp]
# 
# g <- c(up.gene,dn.gene)
# median.value <- apply(expr.data[g,names(purity.vec)],1,median)
# e.g <- names(median.value)[median.value >=1]
# 
# .up.gene <- intersect(e.g,up.gene)
# .dn.gene <- intersect(e.g,dn.gene)

# ############ DOT1L regulated gene ##########
# DOT1L.peak           <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/DOT1L.peak.bed", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
# colnames(DOT1L.peak) <- c('chr','start','end','name')
# DOT1L.peak           <- bedr.sort.region(DOT1L.peak)
# 
# hg19.gene.coordinate           <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/hg19.gene.coordinate.bed", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
# colnames(hg19.gene.coordinate) <- c('chr','strand','start','end','score','gene.name')
# hg19.gene.coordinate           <- hg19.gene.coordinate[,c('chr','start','end','gene.name','score','strand')]
# 
# 
# get.promoter <- function(x) {
#     if(x$strand == '+'){
#         data.frame(chr=x$chr,start = max((x$start) - 1000,0), end = (x$start) + 1000 ,gene.name=x$gene.name) 
#     } else{
#         data.frame(chr=x$chr,start = max((x$end) - 1000,0),   end = (x$end) + 1000, gene.name=x$gene.name) 
#     }   
# }
# hg19.gene.coordinate$idx     <- 1:nrow(hg19.gene.coordinate)
# hg19.gene.promoter           <- ddply(hg19.gene.coordinate,.(idx),get.promoter)
# hg19.gene.promoter           <- hg19.gene.promoter[,c('chr','start','end','gene.name')]
# hg19.gene.promoter$chr       <- as.character(hg19.gene.promoter$chr)
# hg19.gene.promoter$gene.name <- as.character(hg19.gene.promoter$gene.name)
# hg19.gene.promoter           <- unique(hg19.gene.promoter)
# hg19.gene.promoter           <- bedr.sort.region(hg19.gene.promoter,check.chr = FALSE)
# 
# a.int1 <- bedr(input = list(a = hg19.gene.promoter, b = DOT1L.peak), method = "intersect", params = "-loj",check.chr = FALSE)
# DOT1L.regulated.gene <- a.int1$gene.name[a.int1$start.b != -1] %>% as.character() %>% unique


