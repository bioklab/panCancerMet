

up.gene          <- read.csv(file='client-side/output/DE.breast.cancer.R.output/luma.up.csv')$x %>% as.character
require(org.Hs.eg.db)
require(AnnotationDbi)

df <- AnnotationDbi::select(x=org.Hs.eg.db,keys = up.gene,keytype = 'ENSEMBL',columns = c('ENSEMBL','SYMBOL'))
adi.gene <- gene.df$gene.id[gene.df$cell.type == 'Adipocytes']

EMT.gene <- read.table('~/Desktop/EMT.geneset.txt')[,1] %>% as.character()

intersect(df$SYMBOL,EMT.gene)
intersect(df$SYMBOL,adi.gene %>% as.character())

