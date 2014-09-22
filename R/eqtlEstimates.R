setClass("eqtlEstimatesManager",
 representation(fffile="ff_array", call="call", sess="ANY",
        exdate="ANY", shortfac="numeric", geneanno="character", df="numeric",
        summaryList="list"),
        validity=chkeman)

eqtlEstimates = function(smlSet, rhs=~1-1,
   runname="foo", targdir="fooe", 
   geneApply=lapply, 
   shortfac = 10000, checkValid=TRUE, 
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
# store = ff( initdata=0, 
#	dim=c(nsnps, ngenes), 
#	dimnames=list(snpnames, geneNames), 
#  	vmode="short",
#        filename = targff )
# VARIATION: 3 d array
 store = ff( initdata=0, 
	dim=c(nsnps, ngenes, 2), 
	dimnames=list(snpnames, geneNames, c("est.", "s.e.")), 
  	vmode="short",
        filename = targff )
#
 jnk = geneApply( geneNames, function(gene) {
      ex = exprs(smlSet)[gene,]
      fmla = formula(paste("ex", paste(as.character(rhs),collapse=""), collapse=" "))
# VARIATION
      oldwarnopt = options()$warn
      options(warn=-1)
      numans.full = snp.rhs.estimates(fmla, snp.data=snpdata, data=pData(smlSet), 
          family=glmfamily , uncertain=useUncertain)
      options(warn=oldwarnopt)
      numans = sapply(numans.full, "[[", "beta")
# clean if failures
      if (any(badests <- sapply(numans, is.null))) {
                numans[badests] = NA
                numans = unlist(numans)
          }
      numans.var = sapply(numans.full, "[[", "Var.beta")
      if (any(badses <- sapply(numans.var, is.null))) {
                numans.var[badses] = NA
                numans.var = unlist(numans.var)
            }
      store[, gene, 1, add = TRUE] = shortfac * numans
      store[, gene, 2, add = TRUE] = shortfac * sqrt(numans.var)
      NULL
      }) # end gene apply
 close(store)
 exdate = date()
 new("eqtlEstimatesManager", fffile=store, call=theCall, sess=sess, 
        exdate=exdate, shortfac=shortfac, geneanno=annotation(smlSet),
        df=1, summaryList=summfflist)
}

