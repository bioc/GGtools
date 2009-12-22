
setClass("multffManager", contains="list")
setMethod("show", "multffManager", function(object) {
 cat("multffManager instance. The call was:\n")
 print(object$call)
 cat("There are ", length(object$filenames), " ff files.\n")
 cat("Excerpt from first file:\n")
 ngenes = ncol(object$fflist[[1]])
 print(object$fflist[[1]][1:4,1:min(4,ngenes)])
})
 
multffCT = function(listOfSms, gfmla, geneinds=1:10, harmonizeSNPs=FALSE, 
     targdir=".", runname="foo", overwriteFF=TRUE, fillNA=TRUE, ncores=2, vmode="single", ...) {
  theCall = match.call()
  require(ff, quietly=TRUE)
  require(multicore, quietly=TRUE)
  .checkArgsMF( listOfSms, fgmla, geneinds, targdir, runname )
  listOfSms = reduceGenes( listOfSms, geneinds )
  if (harmonizeSNPs) listOfSms = makeCommonSNPs( listOfSms )
  else if(!isTRUE(checkCommonSNPs( listOfSms ))) stop("harmonizeSNPs = FALSE but SNPs not common across listOfSms, run makeCommonSNPs")
  sumScores2ff( listOfSms, gfmla, targdir, runname, theCall, overwriteFF=overwriteFF,
       fillNA=fillNA, write=TRUE, ncores=ncores, vmode=vmode, ... )
}

.checkArgsMF = function( listOfSms, gfmla, geneinds, targdir, runname ) {
  if (!inherits(listOfSms, "list")) stop("listOfSms must inherit from list")
  allc = sapply(listOfSms, function(x) inherits(x, "smlSet"))
  if (!isTRUE(all(allc))) stop("each element of listOfSms must inherit from GGtools smlSet")
 allcsets = sapply(listOfSms, function(x) names(smList(x)))
 if(!is.atomic(allcsets)) stop("probably not all elements of listOfSms have same chromosome set")
 allfn = unique(unlist(lapply(listOfSms, featureNames)))
 allfi = sapply(lapply(listOfSms, featureNames), length)
 if (inherits(geneinds, "character") & !(all(geneinds %in% allfn)))
	stop("some geneinds not in featureNames of listOfSms elements")
 else if (inherits(geneinds, "numeric") & !(max(geneinds)<=max(allfi)))
	stop("some geneinds not in 1:length(featureNames(sms)) for some sms in listOfSms")
 }

reduceGenes = function( listOfSms, geneinds )
  lapply( listOfSms, function(x) x[ geneinds, ] )

makeCommonSNPs = function( listOfSms ) {
  rsidlist = intersectSnps( listOfSms )
  trimSnps( listOfSms, rsidlist )
}

intersectSnps = function( listOfSms ) {
  nsms = length(listOfSms)
  nchr = length(smList(listOfSms[[1]]))
  chnames = names(smList(listOfSms[[1]]))
  rsidlist = lapply(smList(listOfSms[[1]]), colnames)
  for (j in 2:nsms)
    for (k in 1:length(rsidlist))
      rsidlist[[k]]  = intersect(rsidlist[[k]], colnames(smList(listOfSms[[j]])[[k]]))
  rsidlist
}

trimSnps = function( listOfSms, rsidlist ) {
  ans = list()
  cat("making smlSets with common SNPs: \n")
  for (j in 1:length(listOfSms)) {
    cat("smlSet", j, "\nchrom ")
    tmp = listOfSms[[j]]
    for (k in 1:length(rsidlist)) {
      cat(k)
      tmp@smlEnv$smList[[k]] = listOfSms[[j]]@smlEnv$smList[[k]][, rsidlist[[k]]]
      }
    ans[[j]] = tmp
    cat("\n")
    }
  names(ans) = names(listOfSms)
  cat("done.\n")
  gc()
  ans
}

checkCommonSNPs = function( listOfSms ) {
  rsidlist = lapply( smList(listOfSms[[1]]), colnames )
  for (j in 2:length(listOfSms) ) {
    tmp = lapply(smList(listOfSms[[j]]), colnames )
    if (length(tmp) != length(rsidlist)) stop("num chroms not common to all elements of listOfSms")
    for (k in 1:length(rsidlist))
      if (!isTRUE(all.equal(tmp[[k]], rsidlist[[k]]))) {
              cat("sms ", j, " chr ", k, "\n")
              stop("rsid lists not common across smlsets")
              }
     }
  return(TRUE)
}

sumScores2ff = function( listOfSms, gfmla, targdir, runname, theCall=call("1"), 
      overwriteFF=FALSE, fillNA=TRUE, write=TRUE, ncores, vmode, ... ) {
  fnhead = paste(targdir, "/", runname, "_", sep="")
  nsms = length(listOfSms)
  nchr = length(smList(listOfSms[[1]]))
  nsnps = sapply(smList(listOfSms[[1]]), ncol)
  chrnames = names(smList(listOfSms[[1]])) 
  ngenes = length(featn <- featureNames(listOfSms[[1]]))
  generangetag = paste(featn[1], "_",  featn[length(featn)],sep="")
  gnames = featureNames(listOfSms[[1]])
  rslist = lapply(smList(listOfSms[[1]]), colnames)
  filenames = lapply( chrnames, function(x) paste(fnhead, "snpsOnChr", 
      x, "_", generangetag, "_df", nsms, ".ff", sep=""))
  if (write) fflist = lapply(1:nchr, function(x) ff(initdata=0.0, 
    dim=c(nsnps[x], ngenes), 
    dimnames = list(rslist[[x]], gnames),
    vmode=vmode,
    filename=filenames[[x]],
    overwrite=overwriteFF))
  else  fflist = lapply(1:nchr, function(x) ff(initdata=0.0, 
    dim=c(nsnps[x], ngenes), 
    dimnames = list(rslist[[x]], gnames),
    vmode=vmode))
#    filename=paste(fnhead, "snpsOnChr", x, "_", generangetag, ".ff", sep=""),
#    overwrite=overwriteFF)
#  )
  indmat = data.matrix(expand.grid( 1:nchr, 1:nsms ))[,c(2,1)]
  indlist = list()
  for (i in 1:nrow(indmat))
    indlist[[i]] = indmat[i,]
  #for (i in 1:nsms)
  #  mclapply( 1:nchr, function(j) {
  mclapply( indlist, function(indvec) {
    i = indvec[1]
    j = indvec[2]
    #for (j in 1:nchr) {
      cat("sms", i, "chr", j, "\n")
      for (k in 1:ngenes) {
       ex <<- exprs(listOfSms[[i]])[k,]
       gfmla[[2]] = as.name("ex")
       tmpc = snp.rhs.tests( gfmla, snp.data=smList(listOfSms[[i]])[[j]], 
          data=pData(listOfSms[[i]]), family="gaussian", ...)@chisq
       if (fillNA) {
          isna = which(is.na(tmpc))
          if (length(isna)>0) 
             tmpc[isna] = rchisq(length(isna), 1)
          }
       fflist[[j]][,k] = as.ram(fflist[[j]][,k]) + tmpc
       }  # end k
      }, mc.cores=ncores)  # end j/mclapply
   names(fflist) = chrnames
   ans = list(fflist=fflist, call=theCall, runname=runname, targdir=targdir, generangetag=generangetag,
     filenames=filenames, df=nsms)
   assign(runname, new("multffManager", ans))
   save(list=runname, file=paste(runname, ".rda", sep=""))
   invisible(get(runname))
}
  
saveSums = function(sumout, filename="fullout.rda", overwriteFF=TRUE) {
  fflist = sumout$fflist
  chrnames = names(fflist)
  output = list()
  obname = gsub(".rda", "", filename)
  for (i in 1:length(fflist)) {
    obname = paste(sumout$runname, paste("chr", chrnames[i], sep=""), sumout$generangetag, sep="_")
    assign(obname, clone(fflist[[i]], filename=sumout$filenames[[i]], overwrite=overwriteFF))
    output[[i]] = get(obname)
    }
  names(output) = chrnames
  assign(obname, output)
  save(list=obname, file=filename)
}
#
	
# nchroms = length(smList(sms))
# if (inherits(geneinds, "numeric")) genenames = featureNames(sms)[geneinds]
# else genenames = as(geneinds, "character")
# smdims = sapply(smList(sms),dim)
# ngenes = length(geneinds)
# nsamps = length(sampleNames(sms))
# smsname = deparse(substitute(sms))
# sms = sms[ geneinds, ]
# nsnpsPerChr = smdims[2,]
# fns = obnames = rep(NA, nchroms)
# chrnames = names(smList(sms))
# for (i in 1:nchroms) {
#   chunktag = paste(genenames[1], genenames[ngenes], sep="_")
#   basefn = paste(targdir, "/", smsname, ".chr.", chrnames[i], ".", chunktag,
#     sep="")
#   obname = paste(runname, "chr", chrnames[i], chunktag, sep="_")
#   obnames[i] = obname
#   tmp = matrix(rep(0.0, nsnpsPerChr[i]*ngenes), nr=nsnpsPerChr[i])
#   for (j in 1:ngenes) { 
#      ex <<- exprs(sms)[j,]
#      gfmla[[2]] = as.name("ex")
#      ans = snp.rhs.tests( gfmla, snp.data=smList(sms)[[i]], data=pData(sms),
#         family="gaussian", ...)
#      tmp[,j] = ans@chisq
#      gc()
#      }
#   fns[i] = paste(basefn, ".rda", sep="")
#   assign(obname, 
#     ff( initdata=tmp,
#	 overwrite=TRUE,
#         dim=c(nsnpsPerChr[i], ngenes),
#         vmode="double", 
#         dimnames=list( colnames(smList(sms)[[i]]), genenames),
#         filename=paste(basefn, ".ff", sep="")))
#   save(list=obname, file=fns[i])
#   }
#   metadata=data.frame(targdir=targdir, objs=obnames, filenames=fns)
#   list(metadata=metadata, objs=obnames, targdir=targdir, 
#		call=thecall, filenames=fns)
#}
#
#
#  
