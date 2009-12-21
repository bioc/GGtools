
ffcistrans = function(sms, gfmla, geneinds, targdir, runname="foo", ...) {
 require(ff)
 thecall = match.call()
 pvars = names(pData(sms))
 if (inherits(geneinds, "character") & !(all(geneinds %in% featureNames(sms))))
	stop("some geneinds not in featureNames sms")
 else if (inherits(geneinds, "numeric") & !(all(geneinds %in% 1:length(featureNames(sms)))))
	stop("some geneinds not in 1:length(featureNames(sms))")
	
 nchroms = length(smList(sms))
 if (inherits(geneinds, "numeric")) genenames = featureNames(sms)[geneinds]
 else genenames = as(geneinds, "character")
 smdims = sapply(smList(sms),dim)
 ngenes = length(geneinds)
 nsamps = length(sampleNames(sms))
 smsname = deparse(substitute(sms))
 sms = sms[ geneinds, ]
 nsnpsPerChr = smdims[2,]
 fns = obnames = rep(NA, nchroms)
 chrnames = names(smList(sms))
 for (i in 1:nchroms) {
   chunktag = paste(genenames[1], genenames[ngenes], sep="_")
   basefn = paste(targdir, "/", smsname, ".chr.", chrnames[i], ".", chunktag,
     sep="")
   obname = paste(runname, "chr", chrnames[i], chunktag, sep="_")
   obnames[i] = obname
   tmp = matrix(rep(0.0, nsnpsPerChr[i]*ngenes), nr=nsnpsPerChr[i])
   for (j in 1:ngenes) { 
      ex <<- exprs(sms)[j,]
      gfmla[[2]] = as.name("ex")
      ans = snp.rhs.tests( gfmla, snp.data=smList(sms)[[i]], data=pData(sms),
         family="gaussian", ...)
      tmp[,j] = ans@chisq
      gc()
      }
   fns[i] = paste(basefn, ".rda", sep="")
   assign(obname, 
     ff( initdata=tmp,
	 overwrite=TRUE,
         dim=c(nsnpsPerChr[i], ngenes),
         vmode="double", 
         dimnames=list( colnames(smList(sms)[[i]]), genenames),
         filename=paste(basefn, ".ff", sep="")))
   save(list=obname, file=fns[i])
   }
   metadata=data.frame(targdir=targdir, objs=obnames, filenames=fns)
   list(metadata=metadata, objs=obnames, targdir=targdir, 
		call=thecall, filenames=fns)
}
   
#library(GGtools)
#data(hmceuB36.2021)
#ffcistrans = function(sms, gfmla, geneinds, targdir, runname="foo") {
#ut = unix.time(ee <- ffcistrans( hmceuB36.2021, gs~male, 1:50, "."))
#library(multicore)
#ut2 = unix.time( ff <- mclapply( list(1:50, 51:100),
#   function(x) ffcistrans( hmceuB36.2021, gs~male, x, ".")))
 
 
