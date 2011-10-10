setClass("eqtlFDRtab", contains="list")
setMethod("show", "eqtlFDRtab", function(object) {
cat("there were", object$nsnptests, "SNP tested\n", sep=" ")
cat("there were", object$nprobes, "probes\n", sep=" ")
print(object$fdrtab)
cat(object$nc005, "calls at approx FDR = 0.005\n")
cat(object$nc01, "calls at approx FDR = 0.01\n")
cat(object$nc05, "calls at approx FDR = 0.05\n")
#cat(object$nc10, "calls at approx FDR = 0.10\n")
cat("additional elements are:\n", names(object)[-1])
cat("\n")
})

setClass("gwScores", contains="list")
setMethod("show", "gwScores", function(object){
cat("instance of gwScores\n")
cat("elements are:\n", names(object))
cat("\n")
})
setClass("pwScores", contains="list")
setMethod("show", "pwScores", function(object){
cat("instance of pwScores\n")
cat("elements are:\n", names(object))
cat("\n")
})

.policyFDRtab = function(sms, rhs, universe=featureNames(sms),
  policyClo=function(mgr) function(x)topFeats(probeId(x),
     mgr=mgr, ffind=1, n=1), nperm=1, seed=1234, targp=c(.95, .975, .99, .995),
	folderstem="fdrf", geneApply=lapply) {
  set.seed(seed)
  if (nperm > 1) stop("not supporting more than 1 perm at present")
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem)
  per = eqtlTests(permEx(sms), rhs, geneApply=geneApply, 
	targdir=paste("p", folderstem, sep=""))
  policy = policyClo(obs)
  tops = unlist(geneApply(universe, policy))
  policy = policyClo(per)
  ptops = unlist(geneApply(universe, policy))
  nullq = quantile(ptops, targp)
  fcalls = sapply(nullq, function(x)sum(ptops>x))
  scalls = sapply(nullq, function(x)sum(tops>x))
  fdrtab = cbind(pctile=100*targp, thres=nullq, nfalse=fcalls, nsig=scalls, fdr=fcalls/scalls)
  sotops = sort(tops, decreasing=TRUE)
  sptops = sort(ptops, decreasing=TRUE)
  sfdr = sapply(sotops, function(x) sum(sptops>=x)/sum(sotops>=x))
  nf = sfdr*length(sfdr)
  ncall = 1:length(sfdr)
  nc01 = min(which(sfdr >= .01))
  nc05 = min(which(sfdr >= .05))
  nc10 = min(which(sfdr >= .10))
  new("eqtlFDRtab", list(fdrtab=fdrtab, obsmgr=obs, permmgr=per, 
	universe=pm, tops=tops, permtops=ptops,
     	nullq = nullq, targp=targp, ncall=ncall, sfdr=sfdr,
	nc01=nc01, nc05=nc05, nc10=nc10))
}

genewiseScores = function(sms, rhs, targp=c(.95, .975, .99, .995),
	folderstem="fdrf", geneApply=lapply, gene2snpList=NULL) {
#
# factor out the obs and permute steps for genewiseFDRtab, let the 
# permutation occur outside
#
  gn = pm = featureNames(sms)  # gn will be filtered if necessary
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem)
  nsnpsmgd = length(snpsManaged(obs,1))
  nprobesmgd = length(probesManaged(obs,1))
  g2l = gene2snpList  # for revision
  if (is.null(gene2snpList)) {
    tops = unlist(geneApply(pm, function(x)topFeats(probeId(x), mgr=obs, ffind=1, 
		n=1)))   # in this case, tops has names
    }
  else {
    gn = intersect(names(gene2snpList), pm)
    if (length(gn) < 1) stop("gene2snpList is non-null and has null intersection with featureNames(sms)")
    allsn = unlist(lapply(smList(sms), colnames))
    clnsnps = function(x) intersect(x, allsn)  # uses lexical scope to avoid erroneous requests
    g2l = lapply( gene2snpList, clnsnps )
    lens = sapply(g2l, length)
    bad = which(lens==0)  # may lose some genes!
    if (length(bad)>0) {
      g2l = g2l[-bad]
      gn = intersect(names(g2l), pm)
    }
    #tops = sapply(1:length(gn), function(z) {
    #    rs2check = clnsnps(g2l[[ gn[z] ]])
    #    p2check = probeId(gn[z])
    #    max(as.numeric(unlist(obs[rs2check, p2check])))
    #    })
    selector = function(z) {
        rs2check = clnsnps(g2l[[ gn[z] ]])
        p2check = probeId(gn[z])
        tmp = obs[rsid(rs2check), probeId(p2check)]
        nmat = tmp[[1]][,1,drop=FALSE]  # has names, deals with single snp in range
        ind = which.max(as.numeric(unlist(tmp)))
        ans = nmat[ind,1]
        names(ans) = rownames(nmat)[ind]
        ans
        }
    tops = sapply(1:length(gn), selector)
    }
  topdf = data.frame(probes=gn, rsid=names(tops), max.gwscores=tops)
  topdf = topdf[order(topdf$max.gwscores,decreasing=TRUE),]
  scoreq = quantile(tops, targp)
  ndistinctsnps = ifelse(is.null(g2l), nsnpsmgd, length(unique(unlist(g2l))))
  nsnptests = ifelse(is.null(g2l), nsnpsmgd, length(unlist(g2l)))
  nprobes = ifelse(is.null(g2l), nprobesmgd, length(g2l))
  new("gwScores", list(mgr=obs, universe=pm, tops=tops, topdf=topdf,
     scoreq=scoreq, ndistinctsnps=ndistinctsnps, nsnptests = nsnptests,
     nsnpsmgd = nsnpsmgd, nprobesmgd= nprobesmgd,
     nprobes=nprobes))
}

genewiseFDRtab = function(sms, rhs, nperm=1, seed=1234, targp=c(.95, .975, .99, .995),
       folderstem="fdrf", geneApply=lapply, gene2snpList=NULL) {
 thecall = match.call()
 obs = genewiseScores( sms=sms, rhs=rhs, targp=targp, folderstem=folderstem,
	geneApply=geneApply, gene2snpList=gene2snpList )
 set.seed(seed)
 perlist = list()
 for (i in 1:nperm) {
   perlist[[i]] = genewiseScores( sms=permEx(sms), rhs=rhs, targp=targp, 
        folderstem=paste("p_", i, "_",  folderstem, sep=""),
	geneApply=geneApply, gene2snpList=gene2snpList )
   }
# nullq = quantile(per$tops, targp)
 nullq = lapply( perlist, function(x) quantile(x$tops, targp))
# fcalls = sapply(nullq, function(x)sum(per$tops>x))
 fcalls = sapply(1:length(perlist), 
     function(i) sapply(nullq[[i]], function(x)sum(perlist[[i]]$tops>x)))  # need to squelch list character
 fcalls = apply(fcalls, 1, mean)
 nullq = sapply(nullq, function(x)x)
 nullq = apply(nullq,1,mean)
 scalls = sapply(nullq, function(x)sum(obs$tops>x))
##
## do something here to summarize fcalls
##
 fdrtab = cbind(pctile=100*targp, thres=nullq, nfalse=fcalls, nsig=scalls, fdr=fcalls/scalls)
 sotops = sort(obs$tops, decreasing=TRUE)
 sptopslist = list()
 for (i in 1:length(perlist)) {
    sptopslist[[i]] = sort(perlist[[i]]$tops, decreasing=TRUE)
    }
 sptops = apply(sapply(sptopslist, function(x)x), 1, mean)  # check margin here, looks right
 sfdr = sapply(sotops, function(x) sum(sptops>x)/max(c(1,sum(sotops>x))))  # switch to > 10/oct/2011
 nf = sfdr*length(sfdr)
 ncall = 1:length(sfdr)
 nc005 = max(which(sfdr <= .005))
 nc01 = max(which(sfdr <= .01))
 nc05 = max(which(sfdr <= .05))
 nc10 = max(which(sfdr <= .10))
 nc12.5 = max(which(sfdr <= .125))
 nc15 = max(which(sfdr <= .15))
 thresh005 = sptops[nc005]
 thresh01 = sptops[nc01]
 thresh05 = sptops[nc05]
 thresh10 = sptops[nc10]
 thresh12.5 = sptops[nc12.5]
 thresh15 = sptops[nc15]
 tlist = list(thresh005=thresh005, thresh01=thresh01, thresh10=thresh10,
	thresh12.5=thresh12.5,thresh15=thresh15)
 new("eqtlFDRtab", list(fdrtab=fdrtab, obsmgr=obs, permmgr=perlist, 
	unsorted.tops = obs$tops,  topdf=obs$topdf,
	universe=obs$universe, sorted.tops=obs$sotops, sorted.av.permtops=sptops,
        nsnpsmgd = obs$nsnpsmgd, nprobes=obs$nprobes, nsnptests=obs$nsnptests,
     	nullq = nullq, targp=targp, ncall=ncall, sfdr=sfdr,
	nc005=nc005, nc01=nc01, nc05=nc05, nc10=nc10, nc12.5=nc12.5,
        nc15=nc15, threshlist=tlist, thecall=thecall))
}


policywiseFDRtab = function(sms, rhs, nperm=1, seed=1234, targp=c(.95, .975, .99, .995),
       folderstem="fdrf", geneApply=lapply, 
       policyClo=function(mgr) function(x)topFeats(probeId(x),
                mgr=mgr, ffind=1, n=1)) {
 thecall = match.call()
 obs = policywiseScores( sms=sms, rhs=rhs, targp=targp, folderstem=folderstem,
	geneApply=geneApply, policyClo = policyClo )
 set.seed(seed)
 perlist = list()
 for (i in 1:nperm) {
   perlist[[i]] = policywiseScores( sms=permEx(sms), rhs=rhs, targp=targp, 
        folderstem=paste("p_", i, "_",  folderstem, sep=""),
	geneApply=geneApply, policyClo = policyClo )
   }
 nullq = lapply( perlist, function(x) quantile(x$tops, targp))
 fcalls = sapply(1:length(perlist), 
     function(i) sapply(nullq[[i]], function(x)sum(perlist[[i]]$tops>x)))  # need to squelch list character
 fcalls = apply(fcalls, 1, mean)
 nullq = sapply(nullq, function(x)x)
 nullq = apply(nullq,1,mean)
 scalls = sapply(nullq, function(x)sum(obs$tops>x))
##
## do something here to summarize fcalls
##
 fdrtab = cbind(pctile=100*targp, thres=nullq, nfalse=fcalls, nsig=scalls, fdr=fcalls/scalls)
 sotops = sort(obs$tops, decreasing=TRUE)
 sptopslist = list()
 for (i in 1:length(perlist)) {
    sptopslist[[i]] = sort(perlist[[i]]$tops, decreasing=TRUE)
    }
 sptops = apply(sapply(sptopslist, function(x)x), 1, mean)  # check margin here, looks right
 sfdr = sapply(sotops, function(x) sum(sptops>x)/max(c(1,sum(sotops>x))))  # switch to > 10/oct/2011
 nf = sfdr*length(sfdr)
 ncall = 1:length(sfdr)
 nc005 = max(which(sfdr <= .005))
 nc01 = max(which(sfdr <= .01))
 nc05 = max(which(sfdr <= .05))
 nc10 = max(which(sfdr <= .10))
 nc12.5 = max(which(sfdr <= .125))
 nc15 = max(which(sfdr <= .15))
 new("eqtlFDRtab", list(fdrtab=fdrtab, obsmgr=obs, permmgr=perlist, 
	unsorted.tops = obs$tops,  topdf=obs$topdf,
	universe=obs$universe, sorted.tops=obs$sotops, sorted.av.permtops=sptops,
        nsnpsmgd = obs$nsnpsmgd, nprobes=obs$nprobes, nsnptests=obs$nsnptests,
     	nullq = nullq, targp=targp, ncall=ncall, sfdr=sfdr,
	nc005=nc005, nc01=nc01, nc05=nc05, nc10=nc10, nc12.5=nc12.5,
        nc15=nc15, thecall=thecall))
}



policywiseScores = function(sms, rhs, targp=c(.95, .975, .99, .995),
	folderstem="fdrf", geneApply=lapply, 
        policyClo=function(mgr) function(x)topFeats(probeId(x),
                mgr=mgr, ffind=1, n=1)) {
#
# factor out the obs and permute steps for genewiseFDRtab, let the 
# permutation occur outside
#
  gn = pm = featureNames(sms)  # gn will be filtered if necessary
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem)
  policy = policyClo(obs)
  nsnpsmgd = length(snpsManaged(obs,1))
  nprobesmgd = length(probesManaged(obs,1))
  tops = unlist(geneApply(pm, policy))
  topdf = data.frame(probes=gn, rsid=names(tops), max.pwscores=tops)
  topdf = topdf[order(topdf$max.pwscores,decreasing=TRUE),]
  scoreq = quantile(tops, targp)
  g2l = NULL
  ndistinctsnps = ifelse(is.null(g2l), nsnpsmgd, length(unique(unlist(g2l))))
  nsnptests = ifelse(is.null(g2l), nsnpsmgd, length(unlist(g2l)))
  nprobes = ifelse(is.null(g2l), nprobesmgd, length(g2l))
  new("pwScores", list(mgr=obs, universe=pm, tops=tops, topdf=topdf,
     scoreq=scoreq, ndistinctsnps=ndistinctsnps, nsnptests = nsnptests,
     nsnpsmgd = nsnpsmgd, nprobesmgd= nprobesmgd,
     nprobes=nprobes))
}
