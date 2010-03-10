
bestCis = function(ffmgr, slranges, radius=1e6, ffind=1, anno, ncores=10) {
   # get genes in play
# allg = unlist(sapply(ffmgr[[1]], colnames))
 allg = colnames(ffmgr[[1]][[ffind]])
   # get SNPs used in testing
 relevantRS = rownames(ffmgr[[1]][[ffind]])
   # restrict SNPlocs data to relevant snp
 slranges = slranges[ which(slranges$name %in% relevantRS), ]
   # obtain coordinates of genes in play
 gr = geneRanges(allg, anno, extend=radius)
 maxgap <<- 0L
 sspcs = unique(space(slranges))
 if (length(sspcs)>1) warning(paste("slranges included multiple spaces; using", sspcs[1]))
   # find overlaps of gene regions and SNP
 lk = findOverlaps(gr, slranges)[[sspcs[1]]]
 querGnames = allg[ lk@matchMatrix[,1] ]
 indPerGene = split(lk@matchMatrix[,2], querGnames)
   # rs numbers of SNP cis to each gene
 cisrs = mclapply(indPerGene, function(x) slranges$name[x], mc.cores=ncores)
   # get maxchisq of all cis SNP
 names(cisrs) = names(indPerGene)
 wmax = function(x) c(snpind=which.max(x), max=max(x), rsnum=rownames(x)[which.max(x)])
 tmp = mclapply(names(cisrs), function(x) wmax(ffmgr[[1]][[ffind]][cisrs[[x]],x,drop=FALSE]/ffmgr$shortfac))
 tmp = t(matrix(unlist(tmp),nr=3))
 colnames(tmp) = c("snpind", paste("chisq(", ffmgr$df,")", sep=""), "rsnum")
 rownames(tmp) = names(cisrs)  # mclapply does not preserve names
 ans = tmp[match(allg, rownames(tmp) ),]
 ans = data.frame(ans, stringsAsFactors=FALSE, check.names=FALSE)
 ans[,1] = as.numeric(ans[,1])
 ans[,2] = as.numeric(ans[,2])
 ans[,c(3,2)]
}
 
