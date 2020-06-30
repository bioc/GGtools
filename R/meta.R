reduceGenes = function( listOfSms, geneinds )
  lapply( listOfSms, function(x) x[ geneinds, ] )



#eqtlTests = function(smlSet, rhs=~1-1,
#   runname="foo", targdir="foo",
#   geneApply=lapply,
#   shortfac = 100, checkValid=TRUE,
#   useUncertain=TRUE,
#   glmfamily="gaussian") {


meqtlTests = function(listOfSmls, rhslist,
   runname="mfoo", targdir="mfoo", geneApply=lapply, 
   shortfac = 100, checkValid=TRUE, useUncertain=TRUE, 
   glmfamily="gaussian") {
 theCall = match.call()
 sess = sessionInfo()
 geneindex <- 1
 if (missing(glmfamily)) glmfamily="gaussian"
 allfeat = lapply(listOfSmls, featureNames)
 smlSet1 = listOfSmls[[1]]
 fint = allfeat[[1]]
 for (i in 2:length(allfeat)) fint = intersect(fint, allfeat[[i]])
 if (length(fint)==0) stop("null intersection of featureNames for smlSet list elements")
 listOfSmls = reduceGenes( listOfSmls, probeId(fint) )
 listOfSmls = makeCommonSNPs( listOfSmls )

 smlSet1 = listOfSmls[[1]]
 fnhead = paste(targdir, "/", runname, "_", sep="")
 geneNames = featureNames(smlSet1)
 chrName = names(smList(smlSet1))
 cnames = lapply(listOfSmls, function(x) names(smList(x)))
 clens = sapply(cnames,length)
 if (any(clens != 1)) stop("all smlSets must have smList of length 1")
 ngenes = length(geneNames)
 nchr = 1
 if (!file.exists(targdir)) dir.create(targdir)

# there will be one ff file per chromosome which will accumulate
# all information across smlSets

 targff = paste( fnhead, "chr", chrName, ".ff", sep="" )
 allSnpnames = colnames(smList(listOfSmls[[1]])[[1]])
 # do validity check
 jnk = sapply(listOfSmls, validObject)
 if (!all(jnk)) stop("some problem in validity testing after filtering smlSets")
 fffile =
    ff( initdata = 0, dim=c( length(allSnpnames), ngenes),
        dimnames = list(allSnpnames, geneNames), vmode="short",
        filename=targff )
 
# chopped from eqtlTests -- but won't work as such.  hack -- just
# develop summaries on first smlSet in list.  they don't seem to be
# used anyway, except in topFeats for minMAF or minGTF settings...
 summfflist = list()
# if (saveSummaries) {
  # get MAF and minGTF for all SNP
  sumfn = paste(fnhead, chrName, "_summ.ff", sep="")
  if ("multicore" %in% search()) {
    summfflist = geneApply( 1:length(chrName), function(i) ffSnpSummary(smList(smlSet1)[[i]], sumfn[i],
         fac=shortfac))
    } else {
          for (i in 1:length(sumfn))
              summfflist[[chrName[i]]] = ffSnpSummary(smList(smlSet1)[[i]], sumfn[i])
          }
  # ok, now just save references in object
#  }


  store = fffile
  for (theSS in 1:length(listOfSmls)) {
   smlSet = listOfSmls[[theSS]]
   snpdata = smList(smlSet)[[1]]
   gfun = function(gene) {
     geneindex <- geneindex + 1
     ex = exprs(smlSet)[gene,]
     fmla = formula(paste("ex", paste(as.character(rhslist[[theSS]]),collapse=""), collapse=" "))
     numans = snp.rhs.tests(fmla, snp.data=snpdata, 
         data=pData(smlSet), family=glmfamily, uncertain=useUncertain)@chisq
     miss = is.na(numans)
     if (any(miss)) numans[which(miss)] = 0 #rchisq(length(which(miss)), 1)
     store[, gene, add=TRUE] = as.integer(shortfac*numans)
     NULL
    }
    #debug(gfun)
    geneApply( geneNames, gfun ) 
   } # end iterate over smlSet list
  exdate = date()
  new("eqtlTestsManager", fffile=store, call=theCall, sess=sess, 
        exdate=exdate, shortfac=shortfac, geneanno=annotation(smlSet1),
        df=length(listOfSmls), summaryList=summfflist)
}


meta.best.cis.eQTLs.chr = function (smpackvec = c("GGdata", "hmyriB36"), rhslist = list(~1, ~1), folderstem = "mcisScratch",
    radius = 50000, smchr = "20", gchr = "20", schr = "ch20", shortfac=100,
    geneApply = lapply, geneannopk = "illuminaHumanv1.db", snpannopk = snplocsDefault(),
    SMFilterList = list( 
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97),
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97)), 
   exFilterList = list(function(x)x, 
    function(x)x), doPerm=FALSE, excludeRadius=NULL, SSgen=GGBase::getSS
)
{
#
#
    unlink(folderstem, recursive=TRUE)
    cat("get data...")
    smsList = lapply(1:length(smpackvec), function(x) SSgen(
                smpackvec[x], smchr, exFilter=exFilterList[[x]]))
    cat("run smFilter...")   # run earlier than in eqtlTests
    for (i in 1:length(smsList))
         smsList[[i]] =  SMFilterList[[i]](smsList[[i]])
#   get common featurenames
    cat("get common feature names...")
    allpn = lapply(smsList, featureNames)
    commpn = allpn[[1]]
    for (i in 2:length(allpn)) commpn = intersect(commpn, allpn[[i]])
    smsList = lapply(smsList, function(x) x[probeId(commpn),])  # common probes established
    
 
    cat("build map...")
#
# annotation-based list of SNP within radius of coding region
# of gene
#
    cismapObj = getCisMap(radius = radius, gchr = gchr, schr = schr,
        geneannopk = geneannopk, snpannopk = snpannopk, excludeRadius=excludeRadius)
    cismap = namelist(cismapObj)
# map done
# now ensure commonality of mapped probes
#
    allfn = featureNames(smsList[[1]])
    okp = intersect(allfn, names(cismap))
    if (length(okp) < 1)
        stop("no probes common between featureNames and cisMap")
    cat("filter probes in map...")
    cismap = cismap[okp]
    #fsms = sms[probeId(okp),]
    smsList = lapply(smsList, function(x) x[probeId(okp),])
#
# map now agrees with probes in expression component
#
    if (doPerm) {
       cat("permute...")
       smsList = lapply(smsList, permEx)
       }
    cat("tests...")
    mgr = meqtlTests(smsList, rhslist, targdir = folderstem,
        runname = "cis", geneApply = geneApply, shortfac=shortfac)
    mff = fffile(mgr)
    oksn = rownames(mff) 
    cismap = lapply(cismap, function(y)intersect(y,oksn))
    lc = sapply(cismap, length)
    nullc = which(lc == 0)
    if (length(nullc)>0) cismap = cismap[-nullc]
    # genes to use are now names of cismap
    ptested = names(cismap)
    lc = sapply(cismap, length)
    if (length(ptested) == 0) stop("filtering cismap leads to no mapped probes")
    cat("filter...")
    bestcis = geneApply(1:length(ptested), function(pr) {
        curpr = ptested[pr]
        open(mff)
        topind = which.max(mgr[ cismap[[curpr]], curpr])
        bestrs = cismap[[curpr]][topind]
        ans = as.ram(mgr[bestrs, curpr])
        names(ans) = bestrs
        ans
    })
    bestsnp = sapply(bestcis, names)
    names(bestcis) = ptested
    cat("done.\n")
    ans = data.frame(chr=gchr, probe = ptested, snpid = bestsnp, 
        score = as.numeric(bestcis), nsnp = lc,
	stringsAsFactors=FALSE)
    scoredf = ans[order(ans$score, decreasing=TRUE),]
    fullans = GRanges(seqnames=gchr, ranges=cismapObj@generanges[scoredf$probe])
    fullans$score = scoredf$score   # we are assuming that the RangedDat construction does not alter row order!
    fullans$snpid = scoredf$snpid
    fullans$snploc = start(cismapObj@snplocs[scoredf$snpid])
    fullans$radiusUsed = rep(radius, nrow(fullans))
    fullans$nsnp = scoredf$nsnp
    unlink(folderstem, recursive=TRUE)
    fullans
}

meta.best.cis.eQTLs.mchr = function (smpackvec = c("GGdata", "hmyriB36"), rhslist = list(~1,~1), 
      folderstem = "meta.cisScratch", shortfac=100,
    radius = 50000,
    chrnames = as.character(1:22),
    smchrpref = "", gchrpref = "", schrpref = "ch",
    geneApply = lapply,
      geneannopk = "illuminaHumanv1.db",
      snpannopk = snplocsDefault(),
    SMFilterList = list( 
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97),
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97) ),
    exFilterList=list(function(x)x, function(x)x), doPerm=FALSE, excludeRadius=NULL
        ) {
    ans = lapply( chrnames, function(ch) {
            smchr = paste(smchrpref, ch, sep="")
            gchr = paste(gchrpref, ch, sep="")
            schr = paste(schrpref, ch, sep="")
            meta.best.cis.eQTLs.chr(smpackvec = smpackvec, rhslist = rhslist,
             folderstem = folderstem, radius=radius, shortfac=shortfac,
             smchr = smchr, gchr = gchr, schr = schr,
             geneApply = geneApply, geneannopk = geneannopk,
             snpannopk = snpannopk, SMFilterList=SMFilterList,
		exFilterList=exFilterList, doPerm=doPerm, excludeRadius=excludeRadius )
            })
    ans = as(do.call(c, ans), "GRanges")  # GRanges just need c for combination; then mix spaces
    ans[order(elementMetadata(ans)$score, decreasing=TRUE),]
}

meta.best.cis.eQTLs = function(smpackvec = c("GGdata", "hmyriB36"),
   rhslist=list(~1,~1), folderstem="cisScratch", radius=50000, shortfac=100,
    chrnames = as.character(1:22),
    smchrpref = "", gchrpref = "", schrpref = "ch",
    geneApply = lapply,
      geneannopk = "illuminaHumanv1.db",
      snpannopk = snplocsDefault(),
    SMFilterList = list( 
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97),
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97) ), 
  exFilterList=list(function(x)x, function(x)x), nperm=2, excludeRadius=NULL) {
    theCall = match.call()
    obs = meta.best.cis.eQTLs.mchr( smpackvec = smpackvec,
          rhslist=rhslist, folderstem=folderstem, radius=radius, shortfac=shortfac,
          chrnames = chrnames, smchrpref=smchrpref,
          gchrpref=gchrpref, schrpref=schrpref, geneApply=geneApply,
          geneannopk = geneannopk, snpannopk=snpannopk, 
		SMFilterList = SMFilterList, exFilterList=exFilterList, excludeRadius=excludeRadius)
    permans = list()
    alls = numeric()
    fdrs = numeric()
    if (nperm > 0) {
      for (j in 1:nperm) {
        permans[[j]] = meta.best.cis.eQTLs.mchr( smpackvec = smpackvec,
          rhslist=rhslist, folderstem=folderstem, radius=radius, shortfac=shortfac,
          chrnames = chrnames, smchrpref=smchrpref,
          gchrpref=gchrpref, schrpref=schrpref, geneApply=geneApply,
          geneannopk = geneannopk, snpannopk=snpannopk, exFilterList=exFilterList,
          SMFilterList = SMFilterList, doPerm=TRUE, excludeRadius=excludeRadius)  # apply permEx over each filter
        }
      alls = unlist(lapply(permans, function(x)elementMetadata(x)$score))
      fdrs = sapply(elementMetadata(obs)$score, function(x) (sum(alls>=x)/nperm)/sum(elementMetadata(obs)$score>=x))
      elementMetadata(obs)$fdr = fdrs
      obs = obs[order(elementMetadata(obs)$fdr),]
    }  # if nperm == 0, most of following will be empty, but just want scoregr
    new("mcwBestCis", scoregr=obs, allperm=alls, theCall=theCall,
      chromUsed=chrnames, nperm=nperm)
}

