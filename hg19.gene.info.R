hg19.RefSeq.gene.coordinates <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/hg19.RefSeq.gene.coordinates.txt", stringsAsFactors=FALSE)
flag                         <- grepl(x=hg19.RefSeq.gene.coordinates$chrom,pattern = '_')
hg19.RefSeq.gene.coordinates <- hg19.RefSeq.gene.coordinates[!flag,] #Let us remove the genes which are on the un-assembled contigs!
tmp                          <- ddply(hg19.RefSeq.gene.coordinates,.(name2),function(x) data.frame(chrom=x$chrom[1], start=min(x$txStart),end=max(x$txEnd),geneid=x$name2[1],genename=x$name2[1]))
tmp$name2                    <- NULL
hg19.gene.info               <- tmp
hg19.gene.info$chrom         <- gsub(x=hg19.gene.info$chrom,pattern='chr',replacement = '')
hg19.gene.info$chrom         <- as.character(hg19.gene.info$chrom)
hg19.gene.info$geneid        <- as.character(hg19.gene.info$geneid)
hg19.gene.info$genename      <- as.character(hg19.gene.info$genename)
save(file='client-side/output/hg19.gene.info.R.output/hg19.gene.info.RData',list=c('hg19.gene.info'))
