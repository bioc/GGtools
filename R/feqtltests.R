

fsnp.rhs.tests.colch = function(ffref, colch, ...) {
  sdata = new("SnpMatrix", ffref[, colch])
  snp.rhs.tests( ..., snp.data=sdata )
}

concatTests = function (x, y) {
new("GlmTests", snp.names = c(x@snp.names, y@snp.names), chisq = c(x@chisq,
    y@chisq), df = c(x@df, y@df), N = c(x@N, y@N))
}


fsnp.rhs.tests = function(ffref, nchunk=10, applier=lapply, ...) {
  chunks = chunk(1, ncol(ffref), length.out=nchunk)
  ans = applier(chunks, function(ind) fsnp.rhs.tests.colch(ffref, ind, ...))
  if (length(ans)==1) return(unlist(ans))
  res = concatTests( ans[[1]], ans[[2]] )
  if (length(ans)==2) return(res)
  else {
   for (i in 3:length(ans))
      res = concatTests(res, ans[[i]])
   return(res)
  }
}

mfeqtltests = function(listOfEgtSet, listOfRhs=~1, nchunk=10, runname="foo", targdir="foo",
   geneApply=lapply, chromApply=lapply, shortfac=10, scapply=lapply, NA2rchisq=TRUE) {
 thecall = match.call()
 if (length(listOfEgtSet) == 1) stop("intended for use with list of egtSet; call feqtltests directly")
 nset = length(listOfEgtSet)
 ans = feqtltests( listOfEgtSet[[1]], listOfRhs[[1]], nchunk=nchunk, runname=runname, targdir=targdir,
   geneApply=geneApply, chromApply=chromApply, shortfac=shortfac, increment=FALSE, incoming=NULL, scapply=scapply ,
      NA2rchisq=NA2rchisq)
 for (j in 2:length(listOfEgtSet)) {
   tmp = feqtltests( listOfEgtSet[[j]], listOfRhs[[j]], nchunk=nchunk, runname=runname, targdir=targdir,
           geneApply=geneApply, chromApply=chromApply, shortfac=shortfac, increment=TRUE, incoming=ans, scapply=scapply,
           NA2rchisq = NA2rchisq )
   }
 ans@df = nset
 ans@call = thecall
 ans
}
  

feqtltests = function(egtSet, rhs=~1, nchunk=10, runname="foo", targdir="foo",
   geneApply=lapply, chromApply=lapply, shortfac=10, increment=FALSE, incoming=NULL, scapply=lapply, NA2rchisq=TRUE) {
#
# this function works on an out-of-memory representation of genotypes
# there is a double iteration: over chromosomes, and over genes (constant over
#    chromosomes)
# there will be a folder created and, within it, for each chromosome, a single
#    ff file representing chisq scores 
#
 theCall = match.call()
 sess = sessionInfo()
 fnhead = paste(targdir, "/", runname, "_", sep="")
 geneNames = featureNames(egtSet)
 chrNames = names(smList(egtSet))
 if (!file.exists(targdir)) system(paste("mkdir", targdir))
#
# start iteration per chromosome ... set up or receive receptacle
#
 cres = chromApply( chrNames, function(chr) {
      targff = paste( fnhead, "chr", chr, ".ff" , sep="" )
      snpdata = smList(egtSet)[[chr]]
      rownames(snpdata) = sampleNames(egtSet)
      sn = colnames(snpdata) = egtSet@snpNames[[chr]]
      nsnps = ncol(snpdata)
      if (!increment) store = ff( initdata=0, dim=c(nsnps, length(geneNames)), 
                 dimnames=list(sn,geneNames), vmode="short", overwrite=TRUE,
                 filename = targff )
      else store = incoming@fflist[[chr]]
#
# start iteration over genes: populate receptacle store[, gene, add=TRUE] ...
#
      geneApply(geneNames, function(gene) {
           fmla = formula(paste("ex", paste(as.character(rhs),collapse=""), collapse=" "))
           ex = exprs(egtSet)[gene,]
           numans = fsnp.rhs.tests(ffref=snpdata, nchunk=nchunk, applier=scapply,
                        fmla, data=pData(egtSet), 
                        family="Gaussian", uncertain=TRUE )@chisq
#
# to avoid distinguishing NA results (monomorphy in sample) we assign
# from null distribution of test statistics
#
        if (NA2rchisq) {
           if (any(bad <- is.na(numans))) numans[which(bad)] = rchisq(sum(bad), 1)
           }
           store[, gene, add=TRUE] = shortfac*numans
           NULL
           })
      store
      })
  names(cres) = chrNames
  exdate = date()
  new("eqtlTestsManager", fflist=cres, call=theCall, sess=sess, 
        exdate=exdate, shortfac=shortfac, geneanno=annotation(egtSet),
        df=1, summaryList=list())
}

