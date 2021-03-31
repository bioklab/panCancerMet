load('client-side/output/organize.cancer.dependency.R.output/organize.cancer.dependency.RData')

c.cell.line <- intersect(rownames(CCLE_expression.matrix),rownames(Achilles_gene_effect.matrix))
c.gene      <- intersect(colnames(CCLE_expression.matrix),colnames(Achilles_gene_effect.matrix))

idx <- 3
plot(x=CCLE_expression.matrix[c.cell.line[idx],c.gene],y=Achilles_gene_effect.matrix[c.cell.line[idx],c.gene])
lowess(x=CCLE_expression.matrix[c.cell.line[idx],c.gene],y=Achilles_gene_effect.matrix[c.cell.line[idx],c.gene]) %>% lines(lwd=5,col='red')


x <- CCLE_expression.matrix[c.cell.line[idx],c.gene]
y <- Achilles_gene_effect.matrix[c.cell.line[idx],c.gene]
bin.num <- 20
gap     <- max(x) / bin.num

p <- foreach(i = 0:19,.combine='c') %do% {
   flag <- (x >= i*gap) & (x < (i+1)*gap)  
   tmp <- y[flag]
   tmp <- tmp[is.na(tmp) == FALSE]
   sum(tmp < -1) / length(tmp)
   median(tmp)
}
plot(p)




plot(x=CCLE_expression.matrix[c.cell.line,'CDK4'],y=Achilles_gene_effect.matrix[c.cell.line,'CDK4'])
