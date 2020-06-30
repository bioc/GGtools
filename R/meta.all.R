
meta.all.cis.eQTLs.chr = function (minchisq, smpackvec = c("GGdata", "hmyriB36"), rhslist = list(~1, ~1), folderstem = "mcisScratch",
    radius = 50000, smchr = "20", gchr = "20", schr = "ch20", shortfac=100,
    geneApply = lapply, geneannopk = "illuminaHumanv1.db", snpannopk = snplocsDefault(),
    SMFilterList = list( 
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97),
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97)), 
   exFilterList = list(function(x)x, function(x)x), doPerm=FALSE, 
   SSgen=GGBase::getSS
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
        geneannopk = geneannopk, snpannopk = snpannopk)
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
    if (length(ptested) == 0) stop("filtering cismap leads to no mapped probes")

# new below

     ptested = probesManaged(mgr)
     satcis = geneApply(1:length(ptested), function(pr) {
        curpr = ptested[pr]
        oksn = snpsManaged(mgr)
        allscores = as.ram( mgr[, curpr] )
        names(allscores) = oksn
        cisscores = allscores[ intersect(oksn, cismap[[curpr]] ) ]
        if (any(cisscores>=minchisq)) ans = cisscores[ which(cisscores >= minchisq) ]
        else ans = NA
        ans
    })
    bad = unique(c(which(sapply(satcis,length) == 0), which(sapply(satcis,
      function(x)is.na(x[1])))))
    if (length(bad)>0) {
         satcis = satcis[-bad]
         ptested = ptested[-bad] # now a lie, was tested but none survived
         }
    satsnp = lapply(satcis, names)
    satpro = rep(ptested, sapply(satsnp,length))
    cat("done.\n")
    ans = data.frame(chr=gchr, probe = satpro, 
        snpid = unlist(satsnp), score = as.numeric(unlist(satcis)),
        minchisq=minchisq, stringsAsFactors=FALSE)
    scoredf = ans[order(ans$score, decreasing=TRUE),]
# end new

    fullans = GRanges(seqnames=gchr, ranges=cismapObj@generanges[scoredf$probe])
    fullans$score = scoredf$score   # we are assuming that the RangedDat construction does not alter row order!
    fullans$snpid = scoredf$snpid
    fullans$snploc = start(cismapObj@snplocs[scoredf$snpid])
    fullans$radiusUsed = rep(radius, nrow(fullans))
    unlink(folderstem, recursive=TRUE)
    fullans
}

meta.all.cis.eQTLs.mchr = function (minchisq, smpackvec = c("GGdata", "hmyriB36"), rhslist = list(~1,~1), 
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
    exFilterList=list(function(x)x, function(x)x), doPerm=FALSE
        ) {
    ans = lapply( chrnames, function(ch) {
            smchr = paste(smchrpref, ch, sep="")
            gchr = paste(gchrpref, ch, sep="")
            schr = paste(schrpref, ch, sep="")
            meta.all.cis.eQTLs.chr(minchisq, smpackvec = smpackvec, rhslist = rhslist,
             folderstem = folderstem, radius=radius, shortfac=shortfac,
             smchr = smchr, gchr = gchr, schr = schr,
             geneApply = geneApply, geneannopk = geneannopk,
             snpannopk = snpannopk, SMFilterList=SMFilterList,
		exFilterList=exFilterList, doPerm=doPerm )
            })
    ans = as(do.call(c, ans), "GRanges")  # GRanges just need c for combination; then mix spaces
    ans[order(elementMetadata(ans)$score, decreasing=TRUE),]
}

meta.All.cis.eQTLs = function(minchisq, smpackvec = c("GGdata", "hmyriB36"),
   rhslist=list(~1,~1), folderstem="cisScratch", radius=50000, shortfac=100,
    chrnames = as.character(1:22),
    smchrpref = "", gchrpref = "", schrpref = "ch",
    geneApply = lapply,
      geneannopk = "illuminaHumanv1.db",
      snpannopk = snplocsDefault(),
    SMFilterList = list( 
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97),
  function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97) ), 
  exFilterList=list(function(x)x, function(x)x), nperm=2) {
    theCall = match.call()
    obs = meta.all.cis.eQTLs.mchr( minchisq, smpackvec = smpackvec,
          rhslist=rhslist, folderstem=folderstem, radius=radius, shortfac=shortfac,
          chrnames = chrnames, smchrpref=smchrpref,
          gchrpref=gchrpref, schrpref=schrpref, geneApply=geneApply,
          geneannopk = geneannopk, snpannopk=snpannopk, 
		SMFilterList = SMFilterList, exFilterList=exFilterList)
    new("mcwBestCis", scoregr=obs, allperm=numeric(0), theCall=theCall,
      chromUsed=chrnames, nperm=nperm)
}

