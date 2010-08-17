
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
 lk = findOverlaps(gr, slranges)[[as.character(sspcs[1])]]
 querGnames = allg[ lk@matchMatrix[,1] ]
 indPerGene = split(lk@matchMatrix[,2], querGnames)  # some genes will have no overlap
   # rs numbers of SNP cis to each gene
 if (is.loaded("mc_fork", PACKAGE="multicore"))
   cisrs = mclapply(indPerGene, function(x) slranges$name[x], mc.cores=ncores)
 else
   cisrs = lapply(indPerGene, function(x) slranges$name[x] )
   # get maxchisq of all cis SNP
 names(cisrs) = names(indPerGene)
 wmax = function(x) c(snpind=which.max(x), max=max(x), rsnum=rownames(x)[which.max(x)])
 if (is.loaded("mc_fork", PACKAGE="multicore"))
  tmp = mclapply(names(cisrs), function(x) wmax(ffmgr[[1]][[ffind]][cisrs[[x]],x,drop=FALSE]/ffmgr$shortfac))
 else tmp = lapply(names(cisrs), function(x) wmax(ffmgr[[1]][[ffind]][cisrs[[x]],x,drop=FALSE]/ffmgr$shortfac))
 tmp = t(matrix(unlist(tmp),nr=3))
 colnames(tmp) = c("snpind", paste("chisq(", ffmgr$df,")", sep=""), "rsnum")
 rownames(tmp) = names(cisrs)  # mclapply does not preserve names
 ans = tmp[match(names(indPerGene), rownames(tmp) ),]
 ans = data.frame(ans, stringsAsFactors=FALSE, check.names=FALSE)
 ans[,1] = as.numeric(ans[,1])
 ans[,2] = as.numeric(ans[,2])
 ans = ans[,c(3,2)]
 #tmp = list(ans=ans[,c(3,2)], granges=gr, sloc=slranges)
 gstarts = start(gr)
 names(gstarts) = gr$name
 sloc = start(slranges)
 names(sloc) = slranges$name
 snplocs = sloc[ans[,1]]
 gstarts = gstarts[rownames(ans)]
 pv1 = 1-pchisq(ans[,2], ffmgr$df)
 pv2 =pmin(1,2*( 1-pchisq(ans[,2], ffmgr$df)))
 data.frame(gstarts, ans, df=ffmgr$df, snplocs, pv1=pv1, pv2=pv2, check.names=FALSE, stringsAsFactors=FALSE)
}
 
allCisP_1sided = function (ffmgr, slranges, radius = 1e+06, ffind = 1, anno, ncores = 10)
{
    allg = colnames(ffmgr[[1]][[ffind]])
    relevantRS = rownames(ffmgr[[1]][[ffind]])
    slranges = slranges[which(slranges$name %in% relevantRS),
        ]
    gr = geneRanges(allg, anno, extend = radius)
    maxgap <<- 0L
    sspcs = unique(space(slranges))
    if (length(sspcs) > 1)
        warning(paste("slranges included multiple spaces; using",
            sspcs[1]))
    lk = findOverlaps(gr, slranges)[[sspcs[1]]]
    querGnames = allg[lk@matchMatrix[, 1]]
    indPerGene = split(lk@matchMatrix[, 2], querGnames)
    cisrs = mclapply(indPerGene, function(x) slranges$name[x],
        mc.cores = ncores)
    names(cisrs) = names(indPerGene)
    wmax = function(x) c(snpind = which.max(x), max = max(x),
        rsnum = rownames(x)[which.max(x)])
    ans = mclapply(names(cisrs), function(x) 1 - pchisq(ffmgr[[1]][[ffind]][cisrs[[x]],
        x, drop = FALSE]/ffmgr$shortfac, ffmgr$df))
    names(ans) = names(cisrs)
    ans
}

