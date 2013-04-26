chrFilter = function(x, chr="22") {
   anno = annotation(x)
   if (is.null(anno)) stop("no annotation slot")
   require(anno, character.only=TRUE)
   anno = gsub(".db", "", anno)
   chmap = get(paste0(anno, "CHR"))
   print(chr)
   onc = get(chr, revmap(chmap))
   x[ probeId(intersect(featureNames(x), onc)), ]
}

All.cis = 
  function(smpack, rhs=~1, nperm=2, folderstem="cisScratch",
   radius = 50000, shortfac=100, chrnames="22",
   smchrpref="", gchrpref="", schrpref="ch",
   geneApply=lapply, geneannopk = "illuminaHumanv1.db",
   snpannopk = snplocsDefault(), smFilter = function(x)
      nsFilter(MAFfilter(x, lower=.05), var.cutoff=.9), 
   exFilter=function(x)x, keepMapCache=FALSE, SSgen=GGBase::getSS,
   excludeRadius=NULL, estimates=FALSE, ...) {
     thecall = match.call()
     obs = All.cis.mchr( smpack=smpack, rhs=rhs,
      folderstem=folderstem, radius=radius, shortfac=shortfac,
      chrnames=chrnames, smchrpref=smchrpref, gchrpref=gchrpref,
	schrpref=schrpref, geneApply=geneApply, geneannopk=geneannopk,
	snpannopk=snpannopk, smFilter=smFilter, 
	exFilter=exFilter, keepMapCache=keepMapCache, SSgen=SSgen,
        estimates=estimates, ...) 
#
# in following, smFilter is wrapped over permEx
#
     perms = lapply(1:nperm, function(x) All.cis.mchr( smpack=smpack, rhs=rhs,
        folderstem=folderstem, radius=radius, shortfac=shortfac,
        chrnames=chrnames, smchrpref=smchrpref, gchrpref=gchrpref,
	schrpref=schrpref, geneApply=geneApply, geneannopk=geneannopk,
	snpannopk=snpannopk, smFilter=function(x)smFilter(permEx(x)), 
	exFilter=exFilter, keepMapCache=keepMapCache, SSgen=SSgen, estimates=estimates, ...) )
#     list(obs=obs, perms=perms)
     pifdr = function(obs, ps, applier=lapply) {
        nperm=length(ps)/length(obs)
        unlist(applier(obs, function(x) 
           (sum(abs(ps)>=abs(x))/nperm)/sum(abs(obs)>=abs(x)))) }
     obssc = obs$score
     permsc = unlist(lapply(perms, function(x)x$score))
     obs$fdr = pifdr(obssc, permsc)
     obs = obs[order(obs$fdr, -obs$score)]
     smchr = paste0(smchrpref, chrnames)
     obs = bindmaf.simple( smpack, smchr, obs, SSgen, radius )
#bindmaf.simple = function(smpack, smchr, fr, SSgen=GGBase::getSS, rad)
     new("mcwAllCis", obs=obs, perms=perms, theCall=thecall)
}

All.cis.mchr = 
  function(smpack, rhs=~1, folderstem="cisScratch",
   radius = 50000, shortfac=100, chrnames="22",
   smchrpref="", gchrpref="", schrpref="ch",
   geneApply=lapply, geneannopk = "illuminaHumanv1.db",
   snpannopk = snplocsDefault(), smFilter = function(x)
      nsFilter(MAFfilter(x, lower=.05), var.cutoff=.9), 
   exFilter=function(x)x, keepMapCache=FALSE, SSgen=GGBase::getSS,
   excludeRadius=NULL, estimates=FALSE, ...) {
  #
  ans = lapply( chrnames, function(ch) {
    All.cis.chr( smpack=smpack, rhs=rhs,
      folderstem=folderstem, radius=radius, shortfac=shortfac,
      chrname=ch, smchrpref=smchrpref, gchrpref=gchrpref,
	schrpref=schrpref, geneApply=geneApply, geneannopk=geneannopk,
	snpannopk=snpannopk, smFilter=smFilter, 
	exFilter=exFilter, keepMapCache=keepMapCache, SSgen=SSgen, estimates=estimates, ...) } )
  do.call(c, ans)
}
   

All.cis.chr = 
  function(smpack, rhs=~1, folderstem="cisScratch",
   radius = 50000, shortfac=100, chrname="22",
   smchrpref="", gchrpref="", schrpref="ch",
   geneApply=lapply, geneannopk = "illuminaHumanv1.db",
   snpannopk = snplocsDefault(), smFilter = function(x)
      nsFilter(MAFfilter(x, lower=.05), var.cutoff=.9), 
   exFilter=function(x)x, keepMapCache=FALSE, SSgen=GGBase::getSS,
   excludeRadius=NULL, estimates=FALSE, ...) {
#
#
    unlink(folderstem, recursive=TRUE)
    cat("get data...")
    if (length(chrname)>1) stop("chrname must be length 1")
    smchr = paste0(smchrpref, chrname[1])
    sms = SSgen(smpack, smchr, exFilter=exFilter)
    cat("build map...")
#
# annotation-based list of SNP within radius of coding region
# of gene
#
# assumption is that we will compute cis snp to gene mapping ourselves with getCisMap
#
    gchr = paste0(gchrpref, chrname)
    schr = paste0(schrpref, chrname)
    cismapObj = getCisMap(radius = radius, gchr = gchr, schr = schr,
        geneannopk = geneannopk, snpannopk = snpannopk, excludeRadius=excludeRadius)
    cismap = namelist(cismapObj)
#
# prior to may 4 2012 this was run too early... see below
#    if (is.null(mapCache[[gchr]])) mapCache[[gchr]] = cismap  # load cache once for later report, not otherwise reused 4/25/2012
##
## use of gchr here for annotation package
##
#    cat("restrict probes...")
#    if (!is.function(geneannopk)) fsms = restrictProbesToChrom(smFilter(sms), gchr)
### at this point you have cismap which has all relevant probe names
### for chromosome
    cat("run smFilter...")
    sms = smFilter(sms)
    allfn = featureNames(sms)
    okp = intersect(allfn, names(cismap))
    if (length(okp) < 1)
        stop("no probes common between featureNames and cisMap")
    cat("filter probes in map...")
    cismap = cismap[okp]
    fsms = sms[probeId(okp),]
#
# map now agrees with probes in expression component
#
    cat("tests...")
    mgr = eqtlTests(fsms, rhs, targdir = folderstem,
        runname = "cis", geneApply = geneApply, shortfac=shortfac, ...)
    sm = snpsManaged(mgr)
    pl = probesManaged(mgr)
    activeScores = geneApply(pl, function(x)
      mgr[ intersect(sm, cismap[[x]]), x ])
    names(activeScores) = pl
    ntestpg = sapply(activeScores,length)
    gr = cismapObj@generanges[pl]
    gr = rep(gr, ntestpg)
    sn = unlist(lapply(activeScores,rownames))
    gr$snp = sn
    gr$snplocs = start(cismapObj@snplocs[sn])
    gr$score = as.numeric(unlist(activeScores))

    if (estimates) {
       emgr = eqtlEstimates(fsms, rhs, targdir = folderstem,
           runname = "cisEsts", geneApply = geneApply, shortfac=shortfac, ...)
       on.exit(close(emgr@fffile))
       if (!is.open(emgr@fffile)) {
	   open(emgr@fffile)
       }
       edim = dim(emgr@fffile)
       olddim = dim(mgr@fffile)
       if (!all.equal(olddim, edim[1:2])) stop("estimates manager dimension disagrees with test manager dim")
       activeEsts = geneApply(pl, function(x)
         emgr@fffile[ intersect(sm, cismap[[x]]), x, 1 ]/emgr@shortfac)
       activeSEs = geneApply(pl, function(x)
         emgr@fffile[ intersect(sm, cismap[[x]]), x, 2 ]/emgr@shortfac)
       names(activeEsts) = names(activeSEs) = pl
#       ntestpg = sapply(activeEsts,length)
#       gr = cismapObj@generanges[pl]
#       gr = rep(gr, ntestpg)
#       sn = unlist(lapply(activeScores,rownames))
#       gr$snp = sn
#       gr$snplocs = start(cismapObj@snplocs[sn])
       gr$ests = as.numeric(unlist(activeEsts))
       gr$se = as.numeric(unlist(activeSEs))
    }
    #list(mgr=mgr, cismapObj = cismapObj, activeScores=activeScores)
    unlink(folderstem, recursive=TRUE)
    gr
}
     
