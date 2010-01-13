slimCTchr2ff = function (path, chr = 1, dbname = paste("chr", chr, sep = "")) 
{
    require(ff)
    if (!is.numeric(chr)) 
        stop("chr must be numeric index (for Hs, in 1-24)")
    owd = getwd()
    on.exit(setwd(owd))
    tst = try(setwd(path))
    if (inherits(tst, "try-error")) 
        stop("could not change to given folder")
    fn = dir(patt = "\\+.rda$")
    obn = gsub(".rda", "", fn)
    load(fn[1])
    rsn = get(obn[1])$snpnames[[chr]]
    stats = get(obn[1])$chisq[[chr]]
    rm(list = obn[1])
    if (length(obn) > 1) 
        for (i in 2:length(obn)) {
            load(fn[i])
            stats = cbind(stats, get(obn[i])$chisq[[chr]])
            rm(list = obn[i])
        }
    targn = paste(dbname, "chisq",sep="")
    assign(targn, ff(cbind(rsid = rsn, stats), filename=paste(dbname, ".ff", sep=""), dim=c(nrow(stats), ncol(stats)+1),
        dimnames = list(NULL, c("rsid", colnames(stats)))))
    save(list=targn, file=paste(dbname, ".rda", sep=""))
    cat("table written\n")
    NULL
}

ffs2minpBAD = function(ffl, ...) {
 if (!(is(ffl, "list"))) stop("ffl must be a list of ff objects")
 nrs = sapply(ffl, nrow)
 ncs = sapply(ffl, ncol)
 if (!all(ncs == ncs[1])) stop("ffl members must have common ncol, check sapply(ffl, ncol)")
 nff = length(ffl)
 if (nff>1) {
  comrs = ffl[[1]][,1]
  for (i in 2:nff) {
   comrs = intersect(comrs, ffl[[i]][,1])
   }
  clnffl = lapply( ffl, function(x) x[ match(comrs, x[,1]), ] )
  plop1 = function(x){nm=sum(is.na(x));x[is.na(x)] = rchisq(nm,1);x}
  clnffl = lapply( clnffl, plop1)
  clnffl
  summ =  clnffl[[1]]
  for (i in 2:nff) {
    summ = summ + clnffl[[i]]
    }
  }
 else {
  comrs = ffl[[1]][,1]
  summ = ffl[[1]]
  }
 allmax = apply(summ[,-1],1,max)
 # next line, use max for max chisq yield minp for fixed df
 minind = apply(summ[,-1],1, function(x) {tmp = which.max(x); if(length(tmp)==0)tmp=NA;tmp}) # wm can return null
 ans = pmin(1,2*(1-pchisq(allmax,nff)))
 list(rsid = comrs, ans = ans, minind = minind,
        genes = colnames(ffl[[1]])[-1])
}

