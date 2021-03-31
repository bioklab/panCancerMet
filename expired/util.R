### Aim: some functions used by other scripts
require(foreach)
require(dplyr)


#####Function to pick out cell line##########
pick.out.cell.line <- function(expr.of.samples,expr.of.cell.lines,marker.gene){
    marker.gene           <- intersect(rownames(expr.of.samples),(marker.gene))  
    marker.gene           <- intersect(rownames(expr.of.cell.lines),(marker.gene)) 
    correlation.matrix    <- cor(expr.of.samples[marker.gene,],expr.of.cell.lines[marker.gene,],method='spearman')
    cell.line.median.cor  <- apply(correlation.matrix,2,median) %>% sort(decreasing = TRUE)
    best.cell.line        <- names(cell.line.median.cor)[1]
    p.value.vec           <- foreach(cell.line= setdiff(names(cell.line.median.cor),best.cell.line),.combine='c') %do% {
        v                     <- correlation.matrix[,cell.line]
        p.value               <- wilcox.test(correlation.matrix[,best.cell.line],v,alternative = 'greater',paired = TRUE)$p.value
    }
    names(p.value.vec) <- setdiff(names(cell.line.median.cor),best.cell.line)
    fdr.vec            <- p.adjust(p.value.vec,method='fdr')
    list(cell.line.median.cor=cell.line.median.cor,best.cell.line=best.cell.line,compare.fdr.vec=fdr.vec,correlation.matrix = correlation.matrix )
}





