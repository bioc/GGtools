messageCache = new.env(hash=TRUE,parent=emptyenv())

restrictProbesToChrom = function(smlSet, chrom) {
#
# will use CHRLOC for consistency with range computations
#
 ganno = smlSet@annotation
 require(ganno, character.only=TRUE)
 annomappref = gsub(".db", "", ganno)
 map = get(mapn <- paste(annomappref, "CHRLOC", sep=""))
 mk = mappedkeys(map)
 allloc = AnnotationDbi::mget(mk, map, ifnotfound = NA)
 alln = sapply(allloc, function(z) names(z)[1])
 if (any(alln == chrom)) ans = mk[which(alln == chrom)]
 else stop(paste("no CHRLOC elements named with", chrom))
 fn = featureNames(smlSet)
 smlSet[probeId(intersect(fn,ans)),]
}

setMethod("show", "cisMap", function(object) {
   cat("cisMap object with", length(object@namelist), "probes mapped using gene radius", object@radiusUsed, ".\n", sep=" ")
   })
setMethod("namelist", "cisMap", function(cm) cm@namelist)
setMethod("initialize", "cisMap", function(.Object, namelist=list(),
   snplocs=GRanges(), generanges=GRanges(), radiusUsed=numeric()) {
   slot(.Object, "namelist") = namelist
   slot(.Object, "snplocs") = snplocs
   slot(.Object, "generanges") = generanges
   slot(.Object, "radiusUsed") = radiusUsed
   .Object
})
 

.mapSnps2Genes = function(snplocs, generanges, as.GRangesList=FALSE, radiusUsed, ...) {
#
# uses findOverlaps infrastructure to create a list with probes as elements
# and list of snp identifiers as element content
#
   if (is.null(names(generanges))) names(generanges) = as.character(1:length(generanges))
   fo = findOverlaps(generanges, snplocs, ...)
   sn = names(snplocs)
   mm = as.matrix(fo)
   if (prod(dim(mm))==0) stop("no matchMatrix generated in findOverlaps")
#   mmo = mm[order(mm[,2]), ]
   GN = names(generanges)[mm[,1]]
   snpindsByGenes = split(mm[,2], GN) 
#   snpnames = names(snplocs)[mmo[,1]]
# you could get the GRangesList of mapped SNP addresses as
#  lapply( split(snpnames, fac), function(x) snplocs[x] ) 
# here is the mapped list of names
   ans = lapply(snpindsByGenes, function(x)sn[x])
   new("cisMap", namelist=ans, snplocs=snplocs, generanges=generanges,
	radiusUsed=radiusUsed)
}

snpsCisToGenes = function( radius, chr, geneids, genestart, geneend, snpids, snpaddr,
   as.GRangesList=FALSE, ... ) {
  glens = sapply(list(geneids, genestart, geneend), length)
  if (!all(glens==glens[1])) stop("lengths of gene ids, starts, ends must be identical")
  if (!(length(snpids) == length(snpaddr))) stop("lengths of snpids != length snpaddr")
  expranges = IRanges(genestart, geneend)  + radius
  if (any(start(expranges) < 1)) {
    if (!isTRUE(messageCache$negstart)) {
       message(paste("NOTE: expanding gene ranges by radius", radius, "leads to negative start positions that are reset to 1.", sep=" "))
       messageCache$negstart = TRUE
       }
    start(expranges) = pmax(1, start(expranges))
    }
  geneGR = GRanges(seqnames=chr, expranges)
  names(geneGR) = geneids
  snpGR = GRanges(seqnames=chr, IRanges(snpaddr, width=1) )
  names(snpGR) = snpids
  .mapSnps2Genes( snpGR, geneGR, as.GRangesList=as.GRangesList, radiusUsed=radius, ... )
}
  
getCisMap = function( radius=50000, gchr="20", 
  schr="ch20", geneannopk="illuminaHumanv1.db", 
     snpannopk=snplocsDefault(), as.GRangesList=FALSE,
     excludeRadius=NULL ) {
#
# acquires cis map for one chromosome
#
#
# updated 6 Dec 2012: parse geneannopk for FDb in which case we assume
# token has form pkgname:getter
#
  isFDb = FALSE
  if (is.character(geneannopk) && length(grep("^FDb", geneannopk)) == 1) {
    isFDb = TRUE
    toks = strsplit(geneannopk, ":")[[1]]
    geneannopk = toks[1]
    getter = toks[2]
    }
  if (!is.null(excludeRadius) & as.GRangesList) stop("must set as.GRangesList to FALSE when excludeRadius is non-null")
  if (!is.null(excludeRadius) && excludeRadius >= radius) stop("excludeRadius must be < radius")
  if (!is.null(geneannopk) & !is.function(geneannopk)) require(geneannopk, character.only=TRUE, quietly=TRUE)
  require(snpannopk, character.only=TRUE, quietly=TRUE)
#
# 3 jan 2013 -- allow the smpack to provide its own "gene" location
# ranges via featureRanges( gchr ) -- use a function defined in the smpack
#
  if (!isFDb & !is.function(geneannopk)) {
    gpref = gsub(".db", "", geneannopk)
    gcenv = AnnotationDbi::get(paste(gpref, "CHR", sep=""))
    ponc = AnnotationDbi::get(gchr, revmap(gcenv))
    glocenv = AnnotationDbi::get(paste(gpref, "CHRLOC", sep=""))
    glocendenv = AnnotationDbi::get(paste(gpref, "CHRLOCEND", sep=""))
    gstart = AnnotationDbi::mget(ponc, glocenv, ifnotfound=NA)
    gstart = abs(sapply(gstart, "[", 1))
    gend = AnnotationDbi::mget(ponc, glocendenv, ifnotfound=NA)
    gend = abs(sapply(gend, "[", 1))
   } else if (isFDb) {
    featureRanges = get(getter)()
    ponc = names(featureRanges[which(as.character(seqnames(featureRanges)) == gchr)])
    fr = featureRanges[ponc]
    gstart = abs(start(fr))
    gend = abs(end(fr))
   } else if (is.function(geneannopk)) {
      fr = geneannopk(gchr)
      gstart = abs(start(fr))
      gend = abs(end(fr))
      ponc = names(fr)
      } else stop("can't interpret geneannopk")
  bad = is.na(gend) & is.na(gstart)
  if (sum(bad)>0) {
    gstart = gstart[-which(bad)]
    gend = gend[-which(bad)]
    ponc = ponc[-which(bad)]
    }
  slocgr = getSNPlocs(schr, as.GRanges=TRUE)
  if (is.null(names(slocgr))) 
     sids = paste("rs", values(slocgr)$RefSNP_id, sep="")
  else sids = names(slocgr)
  slocs = start(slocgr)
  basic = snpsCisToGenes( radius, gchr, ponc, gstart, gend, sids, slocs, as.GRangesList=as.GRangesList )
  if (is.null(excludeRadius)) return(basic)
  inner = snpsCisToGenes( excludeRadius, gchr, ponc, gstart, gend, sids, slocs, as.GRangesList=as.GRangesList )
  basicn = basic@namelist
  innern = inner@namelist
  basicg = names(basicn)
  innerg = names(innern)
  torii = lapply(1:length(basicg), function(x) {
             curg = basicg[x]
             if (length(innern[[ curg ]])>0)
                         return(setdiff(basicn[[ curg ]], innern[[ curg ]]))
             return(basicn[[curg]])
          })
  names(torii) = basicg
  basic@namelist = torii
  basic@excludeRadius = excludeRadius
  return(basic)
}
 

setMethod("show", "cwBestCis", function(object){
 cat("cwBestCis instance; cited gene ranges are inflated by radius.\n")
 callNextMethod()
})

setMethod("show", "mcwBestCis", function(object) {
 cat("GGtools mcwBestCis instance.  The call was:\n")
 print(object@theCall)
 cat("Best loci for ", np <- length(object@scoregr), " probes are recorded.\n", sep="")
 cat("There were ", object@testCount, " gene:snp tests.\n", sep="")
 toshow = min(np, 4)
 cat("Top ", toshow, "probe:SNP combinations:\n")
 show(object@scoregr[1:toshow])
 cat("====\n")
 cat("use chromsUsed(), fullreport(), etc. for additional information.\n")
})

best.cis.eQTLs.chr = function (smpack = "GGdata", rhs = ~1, folderstem = "cisScratch", shortfac=100,
    radius = 50000, smchr = "20", gchr = "20", schr = "ch20",
    geneApply = lapply, geneannopk = "illuminaHumanv1.db", snpannopk = snplocsDefault(),
    smFilter = function(x) nsFilter(MAFfilter(x, lower = 0.05),
        var.cutoff = 0.97), useME=FALSE, excludeRadius=NULL, exFilter=function(x)x, mapCache=new.env(),
	getDFFITS=FALSE, SSgen=GGBase::getSS)
{
#
#
    unlink(folderstem, recursive=TRUE)
    cat("get data...")
    sms = SSgen(smpack, smchr, exFilter=exFilter)
    cat("build map...")
#
# annotation-based list of SNP within radius of coding region
# of gene
#
# assumption is that we will compute cis snp to gene mapping ourselves with getCisMap
#
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
    if (useME) mgr = eqtlTests.me(fsms, rhs, targdir = folderstem,
        runname = "cis", geneApply = geneApply, shortfac=shortfac)
    else mgr = eqtlTests(fsms, rhs, targdir = folderstem,
        runname = "cis", geneApply = geneApply, shortfac=shortfac)
    mff = fffile(mgr)
    oksn = rownames(mff) 
    cismap = lapply(cismap, function(y)intersect(y,oksn))
    lc = sapply(cismap, length)
    nullc = which(lc == 0)
    if (length(nullc)>0) cismap = cismap[-nullc]
#
# now filter has been applied
#
    if (is.null(mapCache[[gchr]])) mapCache[[gchr]] = cismap  # load cache once for later report, not otherwise reused 4/25/2012
    # genes to use are now names of cismap
    ptested = names(cismap)
    lc = sapply(cismap, length)  # now use this for |s(g)|
    if (length(ptested) == 0) stop("filtering cismap leads to no mapped probes")
    cat("filter...")
    bestcis = lapply(1:length(ptested), function(pr) {
        curpr = ptested[pr]
        topind = which.max(mgr[ cismap[[curpr]], curpr])
        bestrs = cismap[[curpr]][topind]
        ans = as.ram(mgr[bestrs, curpr])
        names(ans) = bestrs
        ans
    })
    bestsnp = sapply(bestcis, names)
    names(bestcis) = ptested
    cat("done.\n")
    ans = data.frame(chr=gchr, probe = ptested, snpid = bestsnp, score = as.numeric(bestcis), nsnp=lc,
	stringsAsFactors=FALSE)
    scoredf = ans[order(ans$score, decreasing=TRUE),]
    fullans = RangedData(seqnames=gchr, ranges=cismapObj@generanges[scoredf$probe])
    fullans$score = scoredf$score   # we are assuming that the RangedDat construction does not alter row order!
    fullans$snpid = scoredf$snpid
    fullans$snploc = start(cismapObj@snplocs[scoredf$snpid])
    fullans$radiusUsed = rep(radius, nrow(fullans))
    fullans$nsnp = scoredf$nsnp
    # dffits option
    if (getDFFITS) fullans$dffits = get.dffits( fsms, rownames(fullans), fullans$snpid )
    unlink(folderstem, recursive=TRUE)
    fullans
}

get.dffits = function( sms, probes, snps ) {
 if (length(probes) != length(snps) ) stop("length probes != length snps")
 numgt = as( smList(sms)[[1]][,snps], "numeric" )
 ex = exprs( sms )[probes,]
 if (!(ncol(numgt) == nrow(ex))) stop("subsetting to probes and snps fails")
 maxdf = function(mod) max(dffits(mod),na.rm=TRUE)
 unlist(sapply(1:nrow(ex), function(i) maxdf(lm(ex[i,,drop=TRUE]~numgt[,i,drop=TRUE]))))
}

best.cis.eQTLs.mchr = function (smpack = "GGdata", rhs = ~1, folderstem = "cisScratch", shortfac=100,
    radius = 50000, 
    chrnames = as.character(1:22),
    smchrpref = "", gchrpref = "", schrpref = "ch",
    geneApply = lapply, 
      geneannopk = "illuminaHumanv1.db", 
      snpannopk = snplocsDefault(),
    smFilter = function(x) nsFilter(MAFfilter(x, lower = 0.05),
        var.cutoff = 0.97), useME=FALSE, 
    excludeRadius=NULL, exFilter=function(x)x, mapCache= new.env(),
    getDFFITS=FALSE, SSgen=GGBase::getSS) {
       ans = lapply( chrnames, function(ch) {
            smchr = paste(smchrpref, ch, sep="")
            gchr = paste(gchrpref, ch, sep="")
            schr = paste(schrpref, ch, sep="")
            best.cis.eQTLs.chr(smpack = smpack, rhs = rhs, 
             folderstem = folderstem, radius=radius, shortfac=shortfac,
             smchr = smchr, gchr = gchr, schr = schr,
             geneApply = geneApply, geneannopk = geneannopk, 
             snpannopk = snpannopk, smFilter=smFilter, exFilter=exFilter,
	     useME=useME, excludeRadius=excludeRadius, mapCache=mapCache,
		getDFFITS=getDFFITS, SSgen=SSgen)
            })
    ans = as(do.call(c, ans), "GRanges")  # RangedData just need c for combination; then mix spaces
    ans[order(elementMetadata(ans)$score, decreasing=TRUE),]
}
    
best.cis.eQTLs = function(smpack = "GGdata",
   rhs=~1, folderstem="cisScratch", radius=50000, shortfac=100,
    chrnames = as.character(1:22),
    smchrpref = "", gchrpref = "", schrpref = "ch",
    geneApply = lapply, 
      geneannopk = "illuminaHumanv1.db", 
      snpannopk = snplocsDefault(),
    smFilter = function(x) nsFilter(MAFfilter(x, lower = 0.05),
        var.cutoff=.97), nperm=2, useME=FALSE, excludeRadius=NULL, exFilter=function(x)x,
	keepMapCache=FALSE, getDFFITS=FALSE, SSgen=GGBase::getSS) {
    theCall = match.call()
    mapCache = new.env()
    obs = best.cis.eQTLs.mchr( smpack = smpack,
          rhs=rhs, folderstem=folderstem, radius=radius, shortfac=shortfac,
          chrnames = chrnames, smchrpref=smchrpref,
	  gchrpref=gchrpref, schrpref=schrpref, geneApply=geneApply,
          geneannopk = geneannopk, snpannopk=snpannopk, smFilter = smFilter, useME=useME, 
		excludeRadius=excludeRadius, exFilter=exFilter, 
		mapCache=mapCache, getDFFITS=getDFFITS, SSgen=SSgen)
    permans = list()
 alls = rep(as.numeric(NA), length(obs))  # initialize in case nperm == 0
  if (nperm > 0) {
    for (j in 1:nperm) {
      permans[[j]] = best.cis.eQTLs.mchr( smpack = smpack,
          rhs=rhs, folderstem=folderstem, radius=radius, shortfac=shortfac,
          chrnames = chrnames, smchrpref=smchrpref,
	  gchrpref=gchrpref, schrpref=schrpref, geneApply=geneApply,
          geneannopk = geneannopk, snpannopk=snpannopk, 
          smFilter = function(x) permEx(smFilter(x)), useME=useME, 
          excludeRadius=excludeRadius, exFilter=exFilter, getDFFITS=getDFFITS,
           SSgen=SSgen)
      }
      alls = unlist(lapply(permans, function(x)elementMetadata(x)$score))
      fdrs = sapply(elementMetadata(obs)$score, function(x) (sum(alls>=x)/nperm)/sum(elementMetadata(obs)$score>=x))
      elementMetadata(obs)$fdr = fdrs
      obs = obs[order(elementMetadata(obs)$fdr),]
    } # end nperm > 0 block
    testCount = length(unlist(as.list(mapCache)))
    if (!keepMapCache) mapCache = new.env()
    new("mcwBestCis", scoregr=obs, allperm=alls, theCall=theCall,
      chromUsed=chrnames, smFilter=smFilter, nperm=nperm, globalMap=mapCache, testCount=testCount)
}

setMethod("chromsUsed", "mcwBestCis", function(x)x@chromUsed)

setMethod("fullreport", c("mcwBestCis", "missing"), function(x, type, ...) {
 x@scoregr
})
setMethod("fullreport", c("mcwBestCis", "character"), function(x, type, ...) {
  if (type == "data.frame") {
   tmp = IRanges::as.data.frame(x@scoregr)
   probeid = rownames(tmp)
   return(data.frame(probeid=probeid, tmp, stringsAsFactors=FALSE))
  } else stop(paste("type", type, "not supported."))
})
  
  
