




xCell.corrected.gene.set <- foreach(x=xCell.data$signatures) %do% {
  x@geneIds
}
names(xCell.corrected.gene.set) <- names(xCell.data$signatures) ############# remove DE genes from xCell signatures. Maybe here we should also remove genes DE between metastatic and primary 

#L <- list()
tmp.df <- foreach(i = names(xCell.corrected.gene.set),.combine='rbind') %do% {
  tmp <- strsplit(x = i,split='%')[[1]] %>% unlist      
  cell.type  <- tmp[1]
  data.source <- tmp[2]
  data.frame(cell.type=cell.type,data.source=data.source,gene.id=xCell.corrected.gene.set[[i]])
  
  #L[[cell.type]] <- c(L[[cell.type]],xCell.corrected.gene.set[[i]])
  #L[[cell.type]] <- unique(L[[cell.type]])
}
#xCell.corrected.gene.set <- L
get.overlapped.gene <- function(x) {
    l <- x$data.source %>% unique %>% length
    ifelse(l >=2,1,0)
}

gene.df <- ddply(tmp.df,.(cell.type,gene.id),get.overlapped.gene)
gene.df <- gene.df[gene.df$V1 == 1,]


