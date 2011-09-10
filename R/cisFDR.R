setClass("eqtlFDRtab", contains="list")
setMethod("show", "eqtlFDRtab", function(object) {
print(object$fdrtab)
cat(object$nc01, "calls at approx FDR = 0.01\n")
cat(object$nc05, "calls at approx FDR = 0.05\n")
cat(object$nc10, "calls at approx FDR = 0.10\n")
cat("additional elements are:\n", names(object)[-1])
cat("\n")
})

.genewiseFDRtab = function(sms, rhs, nperm=1, seed=1234, targp=c(.95, .975, .99, .995),
	folderstem="fdrf", geneApply=mclapply, gene2snpList=NULL) {
  set.seed(seed)
  if (nperm > 1) stop("not supporting more than 1 perm at present")
  pm = featureNames(sms)
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem)
  per = eqtlTests(permEx(sms), rhs, geneApply=geneApply, 
	targdir=paste("p", folderstem, sep=""))
  if (is.null(gene2snpList)) {
    tops = unlist(geneApply(pm, function(x)topFeats(probeId(x), mgr=obs, ffind=1, 
		n=1)))
    ptops = unlist(geneApply(pm, function(x)topFeats(probeId(x), mgr=per, ffind=1, 
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

mgenewiseFDRtab = function(sms, rhs, nperm=2, seed=1234, targp=c(.95, .975, .99, .995),
	folderstem="fdrf", geneApply=mclapply, gene2snpList=NULL) {
  set.seed(seed)
  pm = featureNames(sms)
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem)
  ptopslist = list()
  for (i in 1:nperm) {
    ppref = paste("p", i, "_", sep="")
    pertemp = eqtlTests(permEx(sms), rhs, geneApply=geneApply, 
	targdir=paste(ppref, folderstem, sep=""))
    ptopslist[[i]] = sapply(1:length(gn), function(z) max(as.numeric(unlist(pertemp[rsid(gene2snpList[[gn[z]]]),
        probeId(gn[z]) ]))))
    }
  if (is.null(gene2snpList)) {
    tops = unlist(geneApply(pm, function(x)topFeats(probeId(x), mgr=obs, ffind=1, 
		n=1)))
    ptops = unlist(geneApply(pm, function(x)topFeats(probeId(x), mgr=per, ffind=1, 
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
	folderstem="fdrf", geneApply=mclapply, gene2snpList=NULL) {
#
# factor out the obs and permute steps for genewiseFDRtab, let the 
# permutation occur outside
#
  pm = featureNames(sms)
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem)
  if (is.null(gene2snpList)) {
    tops = unlist(geneApply(pm, function(x)topFeats(probeId(x), mgr=obs, ffind=1, 
		n=1)))
    }
  else {
    gn = intersect(names(gene2snpList), pm)
    if (length(gn) < 1) stop("gene2snpList is non-null and has null intersection with featureNames(sms)")
    tops = sapply(1:length(gn), function(z) max(as.numeric(unlist(obs[rsid(gene2snpList[[gn[z]]]),
        probeId(gn[z]) ]))))
    }
  scoreq = quantile(tops, targp)
  list(mgr=obs, universe=pm, tops=tops, scoreq=scoreq)
}

genewiseFDRtab = function(sms, rhs, seed=1234, targp=c(.95, .975, .99, .995),
       folderstem="fdrf", geneApply=mclapply, gene2snpList=NULL) {
 obs = genewiseScores( sms=sms, rhs=rhs, targp=targp, folderstem=folderstem,
	geneApply=geneApply, gene2snpList=gene2snpList )
 set.seed(seed)
 per = genewiseScores( sms=permEx(sms), rhs=rhs, targp=targp, folderstem=paste("p", folderstem,
        sep=""),
	geneApply=geneApply, gene2snpList=gene2snpList )
 nullq = quantile(per$tops, targp)
 fcalls = sapply(nullq, function(x)sum(per$tops>x))
 scalls = sapply(nullq, function(x)sum(obs$tops>x))
 fdrtab = cbind(pctile=100*targp, thres=nullq, nfalse=fcalls, nsig=scalls, fdr=fcalls/scalls)
 sotops = sort(obs$tops, decreasing=TRUE)
 sptops = sort(per$tops, decreasing=TRUE)
 sfdr = sapply(sotops, function(x) sum(sptops>=x)/sum(sotops>=x))
 nf = sfdr*length(sfdr)
 ncall = 1:length(sfdr)
 nc01 = min(which(sfdr >= .01))
 nc05 = min(which(sfdr >= .05))
 nc10 = min(which(sfdr >= .10))
 new("eqtlFDRtab", list(fdrtab=fdrtab, obsmgr=obs, permmgr=per, 
	universe=obs$pm, tops=obs$tops, permtops=per$tops,
     	nullq = nullq, targp=targp, ncall=ncall, sfdr=sfdr,
	nc01=nc01, nc05=nc05, nc10=nc10))
}
