     pifdr = function(obs, ps, applier=lapply) {
#
# plug-in FDR -- observed score vector and permuted score vector are inputs
#
        rem = length(ps) %% length(obs)
        if (!(rem==0)) warning("permutation vector length not evenly divisible by observed vector length")
        nperm=length(ps)/length(obs)
        unlist(applier(obs, function(x)
           (sum(abs(ps)>=abs(x))/nperm)/sum(abs(obs)>=abs(x)))) 
}


# eqtlTests is being revised to reduce scope creep.  the purpose
# of eqtlTests is to concisely generate all SNP-gene association
# tests, using ff storage and user-selected approaches to
# concurrent computing

ffSnpSummary = function(sm,fn,fac=100) {
#
# this private function will use ff to store information summarizing genotype
# data prior to use with eqtlTests.  The MAF and min GTF are stored
# as short ints, rescaled by fac to allow control of fractional
# precision on these statistics.
#
 dat = col.summary(sm)
 maf = fac*dat[,"MAF"]
 mingtf = fac*apply(dat[,c(5:7)],1,min,na.rm=TRUE)
 FFOVERWRITE = FALSE
 if (.Platform$OS.type == "windows") FFOVERWRITE = TRUE
 if (file.exists(fn)) {
    warning(paste("found existing", fn, "removing..."))
    unlink(fn, recursive=TRUE)
    return(ff( vmode="short", dim=c(length(maf),2),filename=fn, overwrite=FFOVERWRITE,
     dimnames=list(colnames(sm), c("MAF", "mGTF"))))
 }
 ff(initdata=cbind(maf,mingtf), vmode="short", dim=c(length(maf),2),filename=fn,
     overwrite=FFOVERWRITE,
     dimnames=list(colnames(sm), c("MAF", "mGTF")))
}

eqtlTests = function(smlSet, rhs=~1-1,
   runname="foo", targdir="foo", 
   geneApply=lapply, 
   shortfac = 100, checkValid=TRUE, 
   useUncertain=TRUE, 
   glmfamily="gaussian") {
# record call and session information
 theCall = match.call()
 sess = sessionInfo()
# test the input object for validity
 if (checkValid) {
   tmp = validObject(smlSet)
   }
#
# set up receptacle file identifier and top folder
#
 fnhead = paste(targdir, "/", runname, "_", sep="")
 if (!file.exists(targdir)) dir.create(targdir)
#
# acquire feature identifiers, gene and chromosome namdes and count
#
 geneNames = featureNames(smlSet)
 chrNames = names(smList(smlSet))  # will be one or fail
 if (length(chrNames) > 1) stop("multiple chromosome runs no longer supported.")
 ngenes = length(geneNames)
 nchr = length(chrNames)
#
# start genotype summarization -- once used geneApply but this 
# operation may be fast enough in general to avoid this complication
#
 summfflist = list()
  # get MAF and minGTF for all SNP
 sumfn = paste(fnhead, chrNames, "_summ.ff", sep="")
 summfflist = 
     lapply( 1:length(chrNames), 
         function(i) ffSnpSummary(smList(smlSet)[[i]], sumfn[i], 
         fac=shortfac)) 
#
#
#
 fftargs = paste( fnhead, chrNames, ".ff", sep = "" )
 chkf = sapply(fftargs, file.exists)
 if (any(chkf)) stop("some ff target files already exist, please remove or \n \
	 change targdir or runname setting")
 snpdata = smList(smlSet)[[chrNames]]
 snpnames = colnames(snpdata)
 nsnps = ncol(snpdata)
 targff = paste( fnhead, "chr", chrNames, ".ff" , sep="" )
 FFOVERWRITE = FALSE
 if (.Platform$OS.type == "windows") FFOVERWRITE = TRUE
 store = ff( initdata=0, overwrite=FFOVERWRITE,
	dim=c(nsnps, ngenes), 
	dimnames=list(snpnames, geneNames), 
  	vmode="short",
        filename = targff )
 jnk = geneApply( geneNames, function(gene) {
      ex = exprs(smlSet)[gene,]
      fmla = formula(paste("ex", paste(as.character(rhs),collapse=""), collapse=" "))
      numans = snp.rhs.tests(fmla, snp.data=snpdata, data=pData(smlSet), 
          family=glmfamily , uncertain=useUncertain)@chisq
      miss = is.na(numans)
      if (any(miss)) numans[which(miss)] = 0 #rchisq(length(which(miss)), 1)
      store[, gene, add=TRUE] = shortfac*numans
      NULL
      }) # end gene apply
 close(store)
 exdate = date()
 new("eqtlTestsManager", fffile=store, call=theCall, sess=sess, 
        exdate=exdate, shortfac=shortfac, geneanno=annotation(smlSet),
        df=1, summaryList=summfflist)
}

eqtlTests2 = function(smlSet, rhs=~1-1,
   runname="foo", targdir="foo", 
   geneApply=lapply, 
   shortfac = 100, checkValid=TRUE, 
   useUncertain=TRUE, 
   glmfamily="gaussian", inmgr=NULL) {
 if (!is.null(inmgr) & !is(inmgr, "eqtlTestsManager")) stop("inmgr must be eqtlTestsManager instance")
# record call and session information
 theCall = match.call()
 sess = sessionInfo()
# test the input object for validity
 if (checkValid) {
   tmp = validObject(smlSet)
   }
#
# set up receptacle file identifier and top folder
#
 fnhead = paste(targdir, "/", runname, "_", sep="")
 if (!file.exists(targdir)) dir.create(targdir)
#
# acquire feature identifiers, gene and chromosome namdes and count
#
 geneNames = featureNames(smlSet)
 chrNames = names(smList(smlSet))  # will be one or fail
 if (length(chrNames) > 1) stop("multiple chromosome runs no longer supported.")
 ngenes = length(geneNames)
 nchr = length(chrNames)
#
# start genotype summarization -- once used geneApply but this 
# operation may be fast enough in general to avoid this complication
#
 summfflist = list()
  # get MAF and minGTF for all SNP
 sumfn = paste(fnhead, chrNames, "_summ.ff", sep="")
 summfflist = 
     lapply( 1:length(chrNames), 
         function(i) ffSnpSummary(smList(smlSet)[[i]], sumfn[i], 
         fac=shortfac)) 
#
#
#
 fftargs = paste( fnhead, chrNames, ".ff", sep = "" )
 chkf = sapply(fftargs, file.exists)
 if (is.null(inmgr) & any(chkf)) stop("some ff target files already exist, please remove or \n \
	 change targdir or runname setting")
 snpdata = smList(smlSet)[[chrNames]]
 snpnames = colnames(snpdata)
 nsnps = ncol(snpdata)
 targff = paste( fnhead, "chr", chrNames, ".ff" , sep="" )
 FFOVERWRITE = FALSE
 if (.Platform$OS.type == "windows") FFOVERWRITE = TRUE
 if (is.null(inmgr)) {
   store = ff( initdata=0, overwrite=FFOVERWRITE,
	dim=c(nsnps, ngenes), 
	dimnames=list(snpnames, geneNames), 
  	vmode="short",
        filename = targff )
 } else store = inmgr@fffile
 jnk = geneApply( geneNames, function(gene) {
      ex = exprs(smlSet)[gene,]
      fmla = formula(paste("ex", paste(as.character(rhs),collapse=""), collapse=" "))
      numans = snp.rhs.tests(fmla, snp.data=snpdata, data=pData(smlSet), 
          family=glmfamily , uncertain=useUncertain)@chisq
      miss = is.na(numans)
      if (any(miss)) numans[which(miss)] = 0 #rchisq(length(which(miss)), 1)
      store[, gene, add=TRUE] = shortfac*numans
      NULL
      }) # end gene apply
 close(store)
 exdate = date()
 if (!is.null(inmgr)) return(inmgr)
 new("eqtlTestsManager", fffile=store, call=theCall, sess=sess, 
        exdate=exdate, shortfac=shortfac, geneanno=annotation(smlSet),
        df=1, summaryList=summfflist)
}

