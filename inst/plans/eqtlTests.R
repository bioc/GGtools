
eqtlTests = function(smlSet, rhs=~1-1, 
   runname="foo", targdir="foo", geneApply=lapply, chromApply=lapply, 
   shortfac = 100, ... ) {
 fnhead = paste(targdir, "/", runname, "_", sep="")
 geneNames = featureNames(smlSet)
 chrNames = names(smList(smlSet))
 ngenes = length(geneNames)
 nchr = length(chrNames)
 system(paste("mkdir", targdir))
 cres = chromApply( chrNames, function(chr) {
   snpdata = smList(smlSet)[[chr]]
   snpnames = colnames(snpdata)
   nsnps = ncol(snpdata)
   geneApply( geneNames, function(gene) {
     targff = paste( fnhead, "chr", chr, "_", "g", gene, ".ff" , sep="" )
#     store = ff( initdata=0, dim=c(nsnps, 1), dimnames=c(snpnames, gene), vmode="short",
#                 filename = targff)
     ex = exprs(smlSet)[gene,]
     fmla = formula(paste("ex", as.character(rhs), collapse=" "))
     numans = snp.rhs.tests(fmla, snp.data=snpdata, data=pData(smlSet), family="gaussian", ...)@chisq
     miss = is.na(numans)
     if (any(miss)) numans[which(miss)] = rchisq(length(which(miss)), 1)
     ff( initdata=shortfac*numans, dim=c(nsnps, 1), dimnames=list(snpnames, gene), vmode="short",
                 filename = targff )
     }) # end gene apply
  })  # end chr apply
  cres
}
     
     
     
