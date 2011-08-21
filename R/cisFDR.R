setClass("eqtlFDRtab", contains="list")
setMethod("show", "eqtlFDRtab", function(object) {
print(object$fdrtab)
cat("additional elements are:\n", names(object)[-1])
cat("\n")
})

genewiseFDRtab = function(sms, rhs, nperm=1, seed=1234, targp=c(.95, .975, .99, .995),
	folderstem="fdrf", geneApply=mclapply, gene2snpList=NULL) {
  set.seed(seed)
  if (nperm > 1) stop("not supporting more than 1 perm at present")
  pm = featureNames(sms)
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem)
  per = eqtlTests(permEx(sms), rhs, geneApply=geneApply, 
	targdir=paste("p", folderstem, sep=""))
  if (is.null(gene2snpList)) {
    tops = unlist(mclapply(pm, function(x)topFeats(probeId(x), mgr=obs, ffind=1, 
		n=1)))
    ptops = unlist(mclapply(pm, function(x)topFeats(probeId(x), mgr=per, ffind=1, 
		n=1)))
    }
  else {
    gn = intersect(names(gene2snpList), pm)
    if (length(gn) < 1) stop("gene2snpList is non-null and has null intersection with featureNames(sms)")
    tops = sapply(1:length(gn), function(z) max(as.numeric(unlist(obs[rsid(gene2snpList[[gn[z]]]),
        probeId(gn[z]) ]))))
    ptops = sapply(1:length(gn), function(z) max(as.numeric(unlist(per[rsid(gene2snpList[[gn[z]]]),
        probeId(gn[z]) ]))))
    }
  nullq = quantile(ptops, targp)
  fcalls = sapply(nullq, function(x)sum(ptops>x))
  scalls = sapply(nullq, function(x)sum(tops>x))
  fdrtab = cbind(pctile=100*targp, thres=nullq, nfalse=fcalls, nsig=scalls, fdr=fcalls/scalls)
  new("eqtlFDRtab", list(fdrtab=fdrtab, obsmgr=obs, permmgr=per, 
	universe=pm, tops=tops, permtops=ptops,
     	nullq = nullq, targp=targp))
}

policyFDRtab = function(sms, rhs, universe=featureNames(sms),
  policyClo=function(mgr) function(x)topFeats(probeId(x),
     mgr=mgr, ffind=1, n=1), nperm=1, seed=1234, targp=c(.95, .975, .99, .995),
	folderstem="fdrf", geneApply=mclapply) {
  set.seed(seed)
  if (nperm > 1) stop("not supporting more than 1 perm at present")
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem)
  per = eqtlTests(permEx(sms), rhs, geneApply=geneApply, 
	targdir=paste("p", folderstem, sep=""))
  policy = policyClo(obs)
  tops = unlist(mclapply(universe, policy))
  policy = policyClo(per)
  ptops = unlist(mclapply(universe, policy))
  nullq = quantile(ptops, targp)
  fcalls = sapply(nullq, function(x)sum(ptops>x))
  scalls = sapply(nullq, function(x)sum(tops>x))
  fdrtab = cbind(pctile=100*targp, thres=nullq, nfalse=fcalls, nsig=scalls, fdr=fcalls/scalls)
  new("eqtlFDRtab", list(fdrtab=fdrtab, obsmgr=obs, permmgr=per, 
	universe=pm, tops=tops, permtops=ptops,
     	nullq = nullq, targp=targp))
}
