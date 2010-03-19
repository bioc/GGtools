
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
   #targff = paste( fnhead, "chr", chr, "_", "g", gene, ".ff" , sep="" )
   targff = paste( fnhead, "chr", chr, ".ff" , sep="" )
   snpnames = colnames(snpdata)
   nsnps = ncol(snpdata)
   store = ff( initdata=0, dim=c(nsnps, ngenes), dimnames=list(snpnames, geneNames), vmode="short",
                 filename = targff )
   geneApply( geneNames, function(gene) {
     ex = exprs(smlSet)[gene,]
     fmla = formula(paste("ex", paste(as.character(rhs),collapse=""), collapse=" "))
     numans = snp.rhs.tests(fmla, snp.data=snpdata, data=pData(smlSet), family="gaussian", ...)@chisq
     miss = is.na(numans)
     if (any(miss)) numans[which(miss)] = rchisq(length(which(miss)), 1)
     store[, gene, add=TRUE] = shortfac*numans
     NULL
     }) # end gene apply
  store
  })  # end chr apply
  names(cres) = chrNames
  cres
}
     
     
     
