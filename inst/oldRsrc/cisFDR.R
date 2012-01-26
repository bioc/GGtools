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
	folderstem="fdrf", chromApply=lapply, geneApply=lapply, gene2snpList=NULL, ...) {
#
# factor out the obs and permute steps for genewiseFDRtab, let the 
# permutation occur outside
#
#  gn = pm = featureNames(sms)  # gn will be filtered if necessary
# above seems problematic
#
# you are pushing geneExtents and snpRanges through to eqtlTests in ...
#
  obs = eqtlTests(sms, rhs, chromApply=chromApply, geneApply=geneApply, targdir=folderstem, ...)  # could be filtered relative to sms on basis of ...
  nsnpsmgd = length(snpsManaged(obs,1))
  gn <- pm <- probesManaged(obs,1)
  nprobesmgd = length(pm)
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
    g2l = g2l[ gn ]
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
  topdf = data.frame(probes=gn, rsid=names(tops), max.gwscores=tops)  # fragile, assumes names(tops) are rsid
  topdf = topdf[order(topdf$max.gwscores,decreasing=TRUE),]
  scoreq = quantile(tops, targp)
  ndistinctsnps = ifelse(is.null(g2l), nsnpsmgd, length(unique(unlist(g2l))))
  nsnptests = ifelse(is.null(g2l), nsnpsmgd, length(unlist(g2l)))
  nprobes = ifelse(is.null(g2l), nprobesmgd, length(g2l))
  new("gwScores", list(mgr=obs, universe=pm, tops=tops, gn4tops=gn, topdf=topdf,  # everything should work through topdf
     scoreq=scoreq, ndistinctsnps=ndistinctsnps, nsnptests = nsnptests,   # because metadata are there
     nsnpsmgd = nsnpsmgd, nprobesmgd= nprobesmgd,
     nprobes=nprobes))
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
 all.permtops = unlist(lapply(perlist, function(x)x$tops))
 sfdr = sapply(sotops, function(x) sum(sptops>x)/max(c(1,sum(sotops>x))))  # switch to > 10/oct/2011
 sfdr2 = sapply(sotops, function(x) sum(all.permtops>x)/max(c(1,sum(sotops>x))))  # switch to > 10/oct/2011
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
     	nullq = nullq, targp=targp, ncall=ncall, sfdr=sfdr, sfdr2=sfdr2,
	nc005=nc005, nc01=nc01, nc05=nc05, nc10=nc10, nc12.5=nc12.5,
        nc15=nc15, thecall=thecall))
}




genewiseFDRtab = function(sms, rhs, nperm=1, seed=1234, targp=c(.9,.95, .975, .99, .995),
       folderstem="fdrf", chromApply=lapply, geneApply=lapply, gene2snpList=NULL, ...) {
#
# revised 3 nov 2011 with simpler false call enumeration/averaging for FDR
#
 thecall = match.call()
 obs = genewiseScores( sms=sms, rhs=rhs, targp=targp, folderstem=folderstem,
 	chromApply=chromApply,
	geneApply=geneApply, gene2snpList=gene2snpList, ... )
 set.seed(seed)
 perlist = list()
 for (i in 1:nperm) {
   perlist[[i]] = genewiseScores( sms=permEx(sms), rhs=rhs, targp=targp, 
        folderstem=paste("p_", i, "_",  folderstem, sep=""),
	geneApply=geneApply, gene2snpList=gene2snpList, ... )
   }
 nullq = lapply( perlist, function(x) quantile(x$topdf$max.gwscores, targp))
 fcalls = sapply(1:length(perlist), 
     function(i) sapply(nullq[[i]], function(x)sum(perlist[[i]]$topdf$max.gwscores>x)))  # need to squelch list character
 fcalls = apply(fcalls, 1, mean)
 nullq = sapply(nullq, function(x)x)
 nullq = apply(nullq,1,mean)
 scalls = sapply(nullq, function(x)sum(obs$topdf$max.gwscores>x))
##
## do something here to summarize fcalls
##
 fdrtab = cbind(pctile=100*targp, thres=nullq, nfalse=fcalls, nsig=scalls, fdr=fcalls/scalls)
 soprobeids = obs$topdf$probes[ order(obs$topdf$max.gwscores, decreasing=TRUE) ]  # was already ordered!
 sorsid = obs$topdf$rsid[ order(obs$topdf$max.gwscores, decreasing=TRUE) ]
 sotops = sort(obs$topdf$max.gwscores, decreasing=TRUE)
 names(sotops) = soprobeids
 sptopslist = list()
 for (i in 1:length(perlist)) {
    sptopslist[[i]] = perlist[[i]]$topdf$max.gwscores
    }
 #sptops = apply(sapply(sptopslist, function(x)x), 1, mean)  # check margin here, looks right
 sptops = unlist(sptopslist)
 sfdr = sapply(sotops, function(x) (sum(sptops>x)/nperm)/max(c(1,sum(sotops>x))))  # switch to > 10/oct/2011 # use unlist 3 nov 2011
 nf = sfdr*length(sfdr)
 ncall = 1:length(sfdr)
 nc005 = max(which(sfdr <= .005))
 nc01 = max(which(sfdr <= .01))
 nc05 = max(which(sfdr <= .05))
 nc10 = max(which(sfdr <= .10))
 nc12.5 = max(which(sfdr <= .125))
 nc15 = max(which(sfdr <= .15))
 tailps = 1-targp
 amax = function(x,...) ifelse(length(x)==0, NA, max(x, ...))
 n_at_threshs = sapply(tailps, function(x)  amax(which(sfdr <= x)) )
 threshs = sapply(tailps, function(x) sptops[ amax(which(sfdr <= x)) ])
 ncalls = data.frame(ncall=n_at_threshs, fdr=tailps)
 thresh005 = sptops[nc005]
 thresh01 = sptops[nc01]
 thresh05 = sptops[nc05]
 thresh10 = sptops[nc10]
 thresh12.5 = sptops[nc12.5]
 thresh15 = sptops[nc15]
 fullfdrtab = data.frame(probes=soprobeids, # names(sotops), 
     best.QTLid=sorsid, max.gwscores=as.numeric(sotops), fdr=sfdr)
 fullfdrtab = fullfdrtab[ order(fullfdrtab$fdr), ]
 rownames(fullfdrtab) = NULL
 tlist = list(thresh005=thresh005, thresh01=thresh01, thresh10=thresh10,
	thresh12.5=thresh12.5,thresh15=thresh15)

gro = NULL
gr = obs$mgr@geneExtents
if (length(gr)>0) {  # have some loc info
  gro = gr[ as.character(obs$topdf$probes) ]
  values(gro)$chisq = obs$topdf$max.gwscores
  values(gro)$rsid = as.character(obs$topdf$rsid)
  sl = obs$mgr@snpRanges
  slo = sl[ as.character(values(gro)$rsid) ]
  values(gro)$snploc = start(slo)
  values(gro)$score = pmax(0, pmin(16, -log10(sfdr)))
  }
 extras = list(...)

 new("eqtlFDRtab", list(fdrtab=fdrtab, fullfdrtab=fullfdrtab, obsmgr=obs, permmgr=perlist, 
	unsorted.tops = obs$tops,  topdf=obs$topdf, nperm=nperm,
	universe=obs$universe, sorted.tops=obs$sotops, sorted.all.permtops=sptops,
        nsnpsmgd = obs$nsnpsmgd, nprobes=obs$nprobes, nsnptests=obs$nsnptests,
     	nullq = nullq, targp=targp, ncall=ncall, sfdr=sfdr, threshs=threshs,
        n_at_threshs=n_at_threshs, ncalls=ncalls, gro=gro, gene2snpList=gene2snpList,
	rhs=rhs, geneApply=geneApply, chromApply=chromApply, extras=extras,
     thecall=thecall))
}

policywiseScores = function(sms, rhs, targp=c(.95, .975, .99, .995),
	folderstem="fdrf", geneApply=lapply, gene2snpList=NULL, ...) {
#
# factor out the obs and permute steps for genewiseFDRtab, let the 
# permutation occur outside
#
#  gn = pm = featureNames(sms)  # gn will be filtered if necessary
# above seems problematic
#
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem, ...)  # could be filtered relative to sms on basis of ...
  nsnpsmgd = length(snpsManaged(obs,1))
  gn <- pm <- probesManaged(obs,1)
  nprobesmgd = length(pm)
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
    g2l = g2l[ gn ]
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
  topdf = data.frame(probes=gn, rsid=names(tops), max.gwscores=tops)  # fragile, assumes names(tops) are rsid
  topdf = topdf[order(topdf$max.gwscores,decreasing=TRUE),]
  scoreq = quantile(tops, targp)
  ndistinctsnps = ifelse(is.null(g2l), nsnpsmgd, length(unique(unlist(g2l))))
  nsnptests = ifelse(is.null(g2l), nsnpsmgd, length(unlist(g2l)))
  nprobes = ifelse(is.null(g2l), nprobesmgd, length(g2l))
  new("gwScores", list(mgr=obs, universe=pm, tops=tops, gn4tops=gn, topdf=topdf,  # everything should work through topdf
     scoreq=scoreq, ndistinctsnps=ndistinctsnps, nsnptests = nsnptests,   # because metadata are there
     nsnpsmgd = nsnpsmgd, nprobesmgd= nprobesmgd,
     nprobes=nprobes))
}

#ssm(library(GGtools))
policywiseScores = function(sms, rhs, targp=c(.95, .975, .99, .995),
	folderstem="fdrf", geneApply=lapply, policyClo, ...) {
#
# factor out the obs and permute steps for genewiseFDRtab, let the 
# permutation occur outside
#
# policyClo is a function of an eqtlTestsManager, returns a list
#
#
  obs = eqtlTests(sms, rhs, geneApply=geneApply, targdir=folderstem, ...)  # could be 
          # or enhanced filtered relative to sms on basis of ...
  theSnps = snpsManaged(obs, 1)
  theProbes = probesManaged(obs, 1)
  nsnps = length(theSnps)
  nprobes = length(theProbes)
  new("pwScores", policyClo ( obs ))
}

topSameC.policy = function( geneApply ) function( eqtm ) {
    pm = probesManaged( eqtm, 1 )
    tops = unlist(geneApply(pm, function(x)topFeats(probeId(x), mgr=eqtm, ffind=1, 
		n=1)))   # in this case, tops has names
    sn = names(tops)
    list(gene=pm, rsid=sn, score=tops)
}

# request -- if X is a GRanges with seqlengths and X+p leads to an overhang,
# just truncate the offending ranges with a warning

topCis.policy1 = function( snpGR, geneGR, radius, geneApply ) function( eqtm ) {
    pm = probesManaged( eqtm, 1 )
    sm = snpsManaged( eqtm, 1)
# update the ranges to relevant elements
    geneGR = geneGR[ intersect(pm, names(geneGR)) ]
    names(snpGR) = paste("rs", elementMetadata(snpGR)$RefSNP_id, sep="")
    snpGR = snpGR[ intersect(sm, names(snpGR)) ]
# form list of snp names cis to genes
    ovmm = matchMatrix(findOverlaps(snpGR, geneGR))
    if (nrow(ovmm) == 0) stop("no overlaps between radius-extended geneGR and snpGR")
    genev = names(geneGR)[ovmm[,2]]
    snpv = names(snpGR)[ovmm[,1]]
    gene2snpList = split(snpv, genev)
    gn = names(gene2snpList)
 # harvest score resource 
    tops = unlist(geneApply(1:length(gene2snpList), function(x) 
        max(eqtm[ gene2snpList[[x]],gn[x] ][[1]][,1])))
    sn = names(tops)
    list(gene=pm, rsid=sn, score=tops)
}

topCis.policy2 = function( gene2snpList, geneApply ) function( eqtm ) {
    pm = probesManaged( eqtm, 1 )
    sm = snpsManaged( eqtm, 1)
#
# update gene2snpList to relevant elements
#
    gene2snpList = lapply(gene2snpList, function(x)
     intersect(x, sm)) # split(snpv, genev)
    gene2snpList = gene2snpList[ intersect(names(gene2snpList), pm ) ]
    gene2snpList = gene2snpList[ sapply(gene2snpList, length)>0 ]
    gn = names(gene2snpList)
 # harvest score resource 
    tops = unlist(geneApply(1:length(gene2snpList), function(x) 
        sort(eqtm[ rsid(gene2snpList[[x]]), probeId(gn[x]), drop=FALSE ][[1]][,1],
           decreasing=TRUE)[1]))
    sn = names(tops)
    list(gene=pm, rsid=sn, score=tops)
}


#if (!exists("gsl")) gsl = getGene2SnpList( x20fc, "20", "hg19" )
#
#library(parallel)
#topSameC = topSameC.policy( mclapply )
#topCis50kC = topCis.policy2( gsl, mclapply )
#if (!exists("x20")) x20 = getSS("GGdata", "20")
#if (!exists("x20f")) x20f = nsFilter(x20, var.cutoff=.95)
#if (!exists("x20fc")) x20fc = restrictProbesToChrom( x20f, "20" )
#x20fc
#if (!exists("p1"))p1 = policywiseScores( x20fc, ~male, geneApply=mclapply, policyClo =
#  topSameC )
#if (!exists("p2"))p2 = policywiseScores( x20fc, ~male, geneApply=mclapply, policyClo =
#  topCis50kC )

genewiseScores2 = function(sms, rhs, folderstem="fdrf", geneApply=lapply,
    gene2snpList, ... ) {
    ans = policywiseScores( sms, rhs, geneApply, policyClo = topCis.policy2( gene2snpList,
		geneApply ), ... )
    topdf = data.frame(probes=ans$gene, rsid=ans$rsid, max.gwscores=ans$score)
    list(topdf = topdf)
}
	
genewiseFDRtab2 = function(sms, rhs, nperm=1, seed=1234, targp=c(.9,.95, .975, .99, .995),
       folderstem="fdrf", geneApply=lapply, gene2snpList=NULL, ...) {
#
# revised 3 nov 2011 with simpler false call enumeration/averaging for FDR
#
 thecall = match.call()
 obs = genewiseScores2( sms=sms, rhs=rhs, targp=targp, folderstem=folderstem,
	geneApply=geneApply, gene2snpList=gene2snpList, ... )
 set.seed(seed)
 perlist = list()
 for (i in 1:nperm) {
   perlist[[i]] = genewiseScores2( sms=permEx(sms), rhs=rhs, targp=targp, 
        folderstem=paste("p_", i, "_",  folderstem, sep=""),
	geneApply=geneApply, gene2snpList=gene2snpList, ... )
   }
#
# compute plug-in FDR
#
 soprobeids = obs$topdf$probes[ order(obs$topdf$max.gwscores, decreasing=TRUE) ]  # was already ordered!
 sorsid = obs$topdf$rsid[ order(obs$topdf$max.gwscores, decreasing=TRUE) ]
 sotops = sort(obs$topdf$max.gwscores, decreasing=TRUE)
 names(sotops) = soprobeids
 sptopslist = list()
 for (i in 1:length(perlist)) {
    sptopslist[[i]] = perlist[[i]]$topdf$max.gwscores
    }
 sptops = unlist(sptopslist)
 sfdr = sapply(sotops, function(x) (sum(sptops>x)/nperm)/max(c(1,sum(sotops>x))))  # switch to > 10/oct/2011 # use unlist 3 nov 2011
#
# summarize the tail quantiles of permutation distribution
 amax = function(x,...) ifelse(length(x)==0, NA, max(x, ...))
 n_at_threshs = sapply(tailps, function(x)  amax(which(sfdr <= x)) )
 threshs = sapply(tailps, function(x) sptops[ amax(which(sfdr <= x)) ])
 ncalls = data.frame(ncall=n_at_threshs, fdr=tailps)
# now combine these FDR with relevant location information.  when passed on command line
# it gets bound to eqtlTestsManager instance
#
  gr = obs$mgr@geneExtents
  if (is.null(gr)) stop("eqtlTestsManager did not have a geneExtents GRanges")
  sl = obs$mgr@snpRanges
  if (is.null(sl)) stop("eqtlTestsManager did not have a snpRanges GRanges")
  gro = gr[ as.character(obs$topdf$probes) ]
  values(gro)$chisq = obs$topdf$max.gwscores
  values(gro)$rsid = as.character(obs$topdf$rsid)
  sl = obs$mgr@snpRanges
  slo = sl[ as.character(values(gro)$rsid) ]
  values(gro)$snploc = start(slo)
  values(gro)$score = pmax(0, pmin(16, -log10(sfdr)))
  list(obsrd=gro, allpermtops=sptoplist, thecall=theCall, sess=sessionInfo(),
    nperm=nperm, calltab=ncalls)
}

#> getClass(class(dem8f))
#Class "eqtlFDRSummary" [package "GGtools"]
#
#Slots:
#                                                                              
#Name:  allpermtops       obsrd       calls     theCall        sess      genome
#Class:     numeric     GRanges        list        call         ANY   character
#                                                      
#Name:    cisRadius       nperm       nat05       nat01
#Class:     numeric     numeric     numeric     numeric




