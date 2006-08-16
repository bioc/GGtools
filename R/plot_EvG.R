plot_EvG = function(reset, gene, snpid, anno="hgfocus") {
 gn = getpsid(gene, anno)
 y = exprs(reset)[gn,]
 x = snps(reset)[snpid,]
 plot(x,y,ylab=paste("log", gene, "expression"), xlab=paste("minor allele count,",
  snpid), pch=20)
}
