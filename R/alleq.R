fdr = function(x) elementMetadata(x@scoregr)$fdr
minchisqAtFDR = function(x, maxfdr=0.05) min(elementMetadata(x@scoregr)$score[
    fdr(x) <= maxfdr ])

probesWeqtl = function(x, maxfdr=0.05) {
 names(x@scoregr[ fdr(x) <= maxfdr ])
}

All.cis.eQTLs = function (maxfdr = 0.05, inbestcis = NULL, smpack = "GGdata", 
  rhs = ~1, folderstem = "cisScratch", 
  radius = 50000, shortfac=100, chrnames = as.character(1:22), smchrpref = "", 
  gchrpref = "", schrpref = "ch", geneApply = lapply, 
  geneannopk = "illuminaHumanv1.db", 
  snpannopk = snplocsDefault(),
  smFilter4cis = function(x) nsFilter(MAFfilter(clipPCs(x, 1:10),
  lower = 0.05), var.cutoff = 0.85), 
  smFilter4all = function(x) MAFfilter(clipPCs(x, 1:10),
  lower = 0.05),
  nperm = 2, excludeRadius=NULL, exFilter=function(x)x, SSgen=GGBase::getSS) {
 theCall = match.call()
 exdate = date()
cat("PHASE 1: determining cis threshold...\n")
 if (is.null(inbestcis)) {
  btmp = best.cis.eQTLs( smpack=smpack,
    rhs=rhs, folderstem=folderstem, 
    radius=radius, shortfac=shortfac, chrnames=chrnames, smchrpref=smchrpref,
    gchrpref = gchrpref, schrpref = schrpref, geneApply=geneApply,
    geneannopk = geneannopk, snpannopk = snpannopk,
    smFilter=smFilter4cis, nperm=nperm, exFilter=exFilter, SSgen = SSgen )
    }
 else btmp = inbestcis
 kp = probesWeqtl(btmp, maxfdr=maxfdr)
 if (length(kp)<1) stop("no probes meet FDR criterion")
 minchisq = minchisqAtFDR(btmp, maxfdr=maxfdr)
cat("PHASE 2: extracting associations passing cis threshold...\n")
 RDL = list()
 for (i in 1:length(chrnames)) {
    cat("build map...")
#
# annotation-based list of SNP within radius of coding region
#
    gchr = paste(gchrpref, chrnames[i], sep="")
    schr = paste(schrpref, chrnames[i], sep="")
    smchr = paste(smchrpref, chrnames[i], sep="")
    cismapObj = getCisMap(radius = radius, gchr = gchr, schr = schr,
        geneannopk = geneannopk, snpannopk = snpannopk)
    cismap = GGtools:::namelist(cismapObj)
    cokp = intersect(kp, names(cismap))  # mapped probes w/ cis snp on this chr
#
# reduce set of probes
#
     cat("SSgen/filter the probes ...")
     try(unlink(folderstem, recursive=TRUE))
     curss = smFilter4all(SSgen(smpack, smchr, exFilter=exFilter))
     curss = curss[ probeId( intersect( featureNames(curss), cokp ) ), ]
#
#
     cat("test...")
     tmpt = eqtlTests( curss, rhs, 
        targdir=folderstem, runname="all",
	geneApply=geneApply , shortfac=shortfac)
#
# grab hits
#
     ptested = probesManaged(tmpt)
     satcis = geneApply(1:length(ptested), function(pr) {
        curpr = ptested[pr]
        oksn = snpsManaged(tmpt)
        allscores = as.ram( tmpt[, curpr] )
        names(allscores) = oksn
        cisscores = allscores[ intersect(oksn, cismap[[curpr]] ) ]
        ans = cisscores[ cisscores >= minchisq ]
        ans
    })
    satsnp = lapply(satcis, names)
    satpro = rep(ptested, sapply(satsnp,length))
    cat("done.\n")
    ans = data.frame(chr=gchr, probe = satpro, snpid = unlist(satsnp), score = as.numeric(unlist(satcis)),
        maxfdr=maxfdr, stringsAsFactors=FALSE)
    scoredf = ans[order(ans$score, decreasing=TRUE),]
    tmpr = cismapObj@generanges[scoredf$probe]
    svnm = names(tmpr)
    names(tmpr) = NULL
    fullans = RangedData(seqnames=gchr, ranges=tmpr)
    rownames(fullans) = make.names(svnm, unique=TRUE)
    fullans$score = scoredf$score   # we are assuming that the RangedDat construction does not alter row order!
    fullans$snpid = scoredf$snpid
    fullans$probeid = scoredf$probe
    fullans$snploc = start(cismapObj@snplocs[scoredf$snpid])
    fullans$radiusUsed = rep(radius, nrow(fullans))
    fullans$apxfdr = approx(elementMetadata(btmp@scoregr)$score, fdr(btmp), xout=fullans$score)$y
    fullans$maxfdr = rep(maxfdr, nrow(fullans))
    unlink(folderstem, recursive=TRUE)
    RDL[[i]] = fullans
    }
  fulll = do.call(c, RDL)
  new("allSigCis", fulllist=fulll, bestcis=btmp, chromUsed=chrnames,
          theCall = theCall )
}

setMethod("show", "allSigCis", function(object) {
 nr = nrow(object@fulllist)
 cat("All.cis.eQTL output (first", min(c(4,nr)), " of", nr, "rows):\n")
 show(object@fulllist[1:min(4,nr),])
 cat("the call was:\n")
 fd = fdr(object@bestcis)
 cat("there were ", sum(fd < 0.05), " probes with FDR < 0.05.\n")
 cat("note: gene ranges shown are extended by radius (", object@fulllist$radiusUsed[1], ").\n")
 cat("use getCall(), getBest(), getAll() etc. for more info.\n")
})

getBest = function(x) x@bestcis
getCall = function(x) x@theCall
getAll = function(x) x@fulllist

