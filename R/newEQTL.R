#
# in this file we introduce new formal class structures for
# collections of eqtl test runs
#

#
# an eqtlTestsManager can cover a collection of SNP on different
# chromosomes with a single set of genes
# fflist slot holds a list of ff matrices where rows are SNP and columns are
#     genes
# call, sess, exdate geneanno slots are metadata
# shortfac is the scaling factor used to inflate chisq stats so short integer
#     representation has some precision on division by shortfac
# df is d.f. of chisq stat
#
# if em is an eqtlTestsManager instance then em[rsid, probeId] returns
#     a list of chisq statistics properly rescaled
#

#chkeman = function(object){
# eqtlTestsManager validity test
# allgn = lapply(fflist(object), colnames)
# n1 = allgn[[1]]
# chk = sapply(allgn[-1], function(x)all.equal(x,n1))
# if (!all(chk)) return("fflist colnames not common to all elements")
# if (is.null(names(fflist(object)))) return("fflist elements lack names")
# return(TRUE)
#}


snpIdList = function(x) lapply(fflist(x), rownames)

geneIds = function(x) colnames(fflist(x)[[1]])

permuterm = function(l) {
 lens = sapply(l, length)
 eln = rep(names(l), lens)
 dat = unlist(l)
 names(eln) = dat
 eln
}
 
#snpIdMap = function(ids, x) {
#
# use of match seems slow; want to break when find, see below
#
# silx = snpIdList(x)
# m1 = lapply(silx, function(y){ 
#         ans = match(y, ids)
#         names(ans) = y
#         ans
#         })
# m2 = permuterm(m1)
# names(m2) = unlist(silx)
# split(names(m2[ids]), m2[ids])
#}

snpIdMap = function (ids, x) 
{
    chrn = names(x@fflist)
    snpnames = lapply(GGtools:::fflist(x), rownames)
#    cnames = names(snpnames)
#    findchr = function(x) {
#        if (length(x) > 1) 
#            stop("need scalar input")
#        ind = NA
#        for (i in 1:length(snpnames)) {
#            if (x %in% snpnames[[i]]) {
#                ind = i
#                break
#            }
#        }
#        if (is.na(ind)) stop("an rsid was submitted that is not among the snp names in the smlSet for eqtlTests")
#        cnames[ind]
#    }
#    map = sapply(ids, findchr)
#    names(map) = ids
    map = findchr2( ids, snpnames )
    ans = split(names(map), map)  # names of this list will be cardinal numbers
    names(ans) = chrn[ as.numeric(names(ans)) ]
    ans
}

findchr2 = function(snpnvec, listOfSnpn) {
  # determine element in listOfSnpn which matches snpnvec
  # return named vector of element numbers with snpnames as names
  ans = rep(NA, length(snpnvec))
  names(ans) = snpnvec
  for (i in 1:length(listOfSnpn)) {
    mm = match( listOfSnpn[[i]], snpnvec, nomatch=NA )
    if (sum(!is.na(mm))>0) ans[mm] = i
  }
  if (any(is.na(ans))) stop("an rsid was submitted that is not among the snp names in the smlSet for eqtlTests")
  ans
}


ffSnpSummary = function(sm,fn,fac=100) {
 dat = col.summary(sm)
 maf = fac*dat[,"MAF"]
 mingtf = fac*apply(dat[,c(5:7)],1,min,na.rm=TRUE)
 ff(initdata=cbind(maf,mingtf), vmode="short", dim=c(length(maf),2),filename=fn,
     dimnames=list(colnames(sm), c("MAF", "mGTF")))
}


 
eqtlTests = function(smlSet, rhs=~1-1,
   runname="foo", targdir="foo", geneApply=lapply, chromApply=lapply,
   shortfac = 100, computeZ=FALSE, checkValid=TRUE, saveSummaries=TRUE, uncert=TRUE, family, genegran=50, ... ) {
 theCall = match.call()
 if (checkValid) {
   tmp = validObject(smlSet)
   }
 if (missing(family)) family="gaussian"
 geneindex <<- 1
 sess = sessionInfo()
 fnhead = paste(targdir, "/", runname, "_", sep="")
 geneNames = featureNames(smlSet)
 chrNames = names(smList(smlSet))
 ngenes = length(geneNames)
 nchr = length(chrNames)
 if (!file.exists(targdir)) system(paste("mkdir", targdir))
 summfflist = list()
 if (saveSummaries) {
  # get MAF and minGTF for all SNP
  sumfn = paste(fnhead, chrNames, "_summ.ff", sep="")
  if ("multicore" %in% search()) {
    summfflist = mclapply( 1:length(chrNames), function(i) ffSnpSummary(smList(smlSet)[[i]], sumfn[i], 
         fac=shortfac)) 
    } else {
          for (i in 1:length(sumfn))
              summfflist[[chrNames[i]]] = ffSnpSummary(smList(smlSet)[[i]], sumfn[i])
          }
  # ok, now just save references in object
  }
 cres = chromApply( chrNames, function(chr) {
   snpdata = smList(smlSet)[[chr]]
   #targff = paste( fnhead, "chr", chr, "_", "g", gene, ".ff" , sep="" )
   targff = paste( fnhead, "chr", chr, ".ff" , sep="" )
   snpnames = colnames(snpdata)
   nsnps = ncol(snpdata)
   store = ff( initdata=0, dim=c(nsnps, ngenes), dimnames=list(snpnames, geneNames), vmode="short",
                 filename = targff )
   geneApply( geneNames, function(gene) {
     if (options()$verbose & geneindex %% genegran == 0) cat(gene, "..")
     geneindex <<- geneindex + 1
     if (options()$verbose & geneindex %% 8*genegran == 0) cat("\n")
     ex = exprs(smlSet)[gene,]
     fmla = formula(paste("ex", paste(as.character(rhs),collapse=""), collapse=" "))
     numans = snp.rhs.tests(fmla, snp.data=snpdata, data=pData(smlSet), 
         family=family , ...)@chisq
#uncertain=uncert
     if (computeZ) {
       numans = sqrt(numans)
       signl = snp.rhs.estimates( fmla, snp.data=snpdata, data=pData(smlSet), family=family, ... )
       bad = which(unlist(lapply(signl, is.null)))
       if (length(bad)>0) signl[bad] = list(beta=NA)
       ifelse(unlist(signl)>=0, 1, -1)
       numans = numans*signl
     }
     miss = is.na(numans)
     if (any(miss) & !computeZ) numans[which(miss)] = rchisq(length(which(miss)), 1)
     if (any(miss) & computeZ) numans[which(miss)] = rnorm(length(which(miss)))
     store[, gene, add=TRUE] = shortfac*numans
     NULL
     }) # end gene apply
  store
  })  # end chr apply
  names(cres) = chrNames
  exdate = date()
  new("eqtlTestsManager", fflist=cres, call=theCall, sess=sess, 
        exdate=exdate, shortfac=shortfac, geneanno=annotation(smlSet),
        df=1, summaryList=summfflist)
}

# director for group of managers
#
#chkmgrs = function(object) {
#   mcl = sapply(mgrs(object), class)
#   chkc = sapply(mgrs(object), function(x) is(x, "eqtlTestsManager"))
#   if (!all(chkc)) return("mgrs slot must only contain list of entities inheriting from eqtlTestsManager")
#   annos = sapply(mgrs(object), function(x)x@geneanno)
#   if (!all(annos==annos[1])) return("managers do not have identical gene annotation source")
#   sids = lapply(mgrs(object), snpIdList)
#   slchk = sapply(sids, function(x) all.equal(x, sids[[1]]))
#   if (!all(sapply(slchk,isTRUE))) return("managers do not have identical SNP lists")
#   return(TRUE)
#}
#   

mkCisTransDirector = function(dl, indexdbname, snptabname, probetabname, probeanno, commonSNPs=TRUE) {
 cd = new("cisTransDirector", indexdbname=indexdbname, shortfac=shortfac(dl[[1]]), mgrs=dl,
     snptabname=snptabname, probetabname=probetabname, probeanno=probeanno)
 ffrefs = mkDirectorDb(cd, commonSNPs)
 cd@snptabref = ffrefs$snptabref
 cd@probetabref = ffrefs$probetabref
 cd
}

mkDirectorDb = function(cd, commonSNPs=TRUE) {
#
# objective here is a small footprint dump to two ff files that
# will serve as indexes
# [indexdbname]_snpnames.ff will map from snpids to chr
# [indexdbname]_probenames.ff will have all gene names and the manager indices
#

 vecs2ff = function(nmdlist, filename) {
# support for dumping index data
   vlist = names(nmdlist)
   ref = ff(initdata=nmdlist[[2]], file=filename, dim=c(length(nmdlist[[1]]), length(nmdlist)-1),
     vmode="short")
   rownames(ref) = as.character(nmdlist[[1]])
   ref
  }


#
#
 if (commonSNPs) {
   f1 = fflist(mgrs(cd)[[1]])
   rsids = unlist(lapply(fflist(mgrs(cd)[[1]]), rownames))
# cn here denotes chromosome names.  but we want our
# ff to be populated with short ints as indices.  so we
# will use integers ... should work with fflist indexing
   cn = 1:length(names(f1)) #names(f1)
   chrs = rep(cn, sapply(f1, nrow))
   mgr = rep(1, length(rsids))
   snptabref = vecs2ff( list(snpid=rsids, chr=chrs),  
       paste(cd@indexdbname, "_", cd@snptabname, ".ff", sep="") )
 }
 else stop("only handling managers with common SNP fields")

 allg = lapply(mgrs(cd), function(x) colnames(fflist(x)[[1]])) # 
 pids = allg[[1]]
 mgr = rep(1, length(pids))
 if (length(allg) > 1) for (i in 2:length(allg)) {
                           pids = c(pids, allg[[i]])
                           mgr = c(mgr, rep(i, length(allg[[i]])))
                       }
  
  probetabref = vecs2ff( list(probeid=pids, mgr=mgr), 
            paste(cd@indexdbname, "_", cd@probetabname, ".ff", sep=""))
  list(snptabref=snptabref, probetabref=probetabref)
}


ieqtlTests = function (smlSet, rhs = ~1 - 1, rules, runname = "ifoo", targdir = "ifoo", 
    geneApply = lapply, chromApply = lapply, shortfac = 100, 
    computeZ = FALSE, uncert=TRUE, saveSummaries=TRUE,
    family, ...) 
{
    theCall = match.call()
    if (missing(family)) family="gaussian"
    sess = sessionInfo()
    fnhead = paste(targdir, "/", runname, "_", sep = "")
    geneNames = featureNames(smlSet)
    chrNames = names(smList(smlSet))
    ngenes = length(geneNames)
    nchr = length(chrNames)
    system(paste("mkdir", targdir))
#
# following just grabbed from eqtlTests
#
 summfflist = list()
if (saveSummaries) {
  # get MAF and minGTF for all SNP
  sumfn = paste(fnhead, chrNames, "_summ.ff", sep="")
  if ("multicore" %in% search()) {
    summfflist = mclapply( 1:length(chrNames), function(i) ffSnpSummary(smList(smlSet)[[i]], sumfn[i],
         fac=shortfac))
    } else {
          for (i in 1:length(sumfn))
              summfflist[[chrNames[i]]] = ffSnpSummary(smList(smlSet)[[i]], sumfn[i])
          }
  # ok, now just save references in object
  }

    cres = chromApply(chrNames, function(chr) {
        snpdata = smList(smlSet)[[chr]]
        targff = paste(fnhead, "chr", chr, ".ff", sep = "")
        snpnames = c(colnames(snpdata), names(rules))
        nsnps = length(snpnames)
        store = ff(initdata = 0, dim = c(nsnps, ngenes), dimnames = list(snpnames, 
            geneNames), vmode = "short", filename = targff)
        geneApply(geneNames, function(gene) {
            ex = exprs(smlSet)[gene, ]
            fmla = formula(paste("ex", paste(as.character(rhs), 
                collapse = ""), collapse = " "))
            numans = snp.rhs.tests(fmla, snp.data = snpdata, 
                data = pData(smlSet), family = family, uncertain=uncert, 
                ...)@chisq
            numansi = snp.rhs.tests(fmla, snp.data = snpdata, uncertain=uncert,
                data = pData(smlSet), family = family, rules = rules, 
                ...)@chisq
            numans = c(numans, numansi)
            if (computeZ) {
                numans = sqrt(numans)
                signl = snp.rhs.estimates(fmla, snp.data = snpdata, 
                  data = pData(smlSet), family = family, 
                  ...)
                bad = which(unlist(lapply(signl, is.null)))
                if (length(bad) > 0) 
                  signl[bad] = list(beta = NA)
                ifelse(unlist(signl) >= 0, 1, -1)
                numans = numans * signl
            }
            miss = is.na(numans)
            if (any(miss) & !computeZ) 
                numans[which(miss)] = rchisq(length(which(miss)), 
                  1)
            if (any(miss) & computeZ) 
                numans[which(miss)] = rnorm(length(which(miss)))
            store[, gene, add = TRUE] = shortfac * numans
            NULL
        })
        store
    })
    names(cres) = chrNames
    exdate = date()
    new("eqtlTestsManager", fflist = cres, call = theCall, sess = sess, 
        exdate = exdate, shortfac = shortfac, geneanno = annotation(smlSet), 
        df=1, summaryList=summfflist)
}

getNamedLocs = function(slpack="SNPlocs.Hsapiens.dbSNP.20100427", chrtok) {
 require(slpack, character.only=TRUE)
 if (slpack == "SNPlocs.Hsapiens.dbSNP.20100427" && length(grep("chr", chrtok))>0) {
    chrtok = gsub("chr", "ch", chrtok)
    warning("don't use chrNN with 20100427 snplocs package ... trying chNN ...")
    }
 locdf = getSNPlocs(chrtok)
 rsid = paste("rs", locdf$RefSNP_id, sep="")
 locs = locdf$loc
 names(locs) = rsid
 locs
}
 

getGRanges = function(mgr, ffind, geneind, seqnames, namedlocs) {
  if (length(geneind) != 1) stop("geneind must be scalar")
  snpids = snpIdList(mgr)[[ffind]]
  scores = fflist(mgr)[[ffind]][, geneind]/shortfac(mgr)
  names(scores) = snpids
  okids = intersect(names(namedlocs), snpids)
  oklocs = namedlocs[okids]
  okscores = scores[okids]
  n = length(okids)
  tmp = GRanges(IRanges(oklocs, width=1), seqnames=rep(seqnames,n), 
        score=as.numeric(-log10(1-pchisq(okscores, mgr@df))),
        chisq=as.numeric(okscores), df=rep(mgr@df, n))
  names(tmp) = okids
  tmp
}

cisRanges = function(probeids, chr, anno, radius=5e5, useEnd=FALSE) {
 require(GenomicRanges)
 require(anno, character.only=TRUE)
 goodchr = gsub("chr", "", chr)
 chrs = mget( probeids, get(paste(gsub(".db", "", anno), "CHR", sep="")), ifnotfound=NA)
 thechrs = unique(unlist(na.omit(chrs)))
 if (length(thechrs)>1) stop("probeids supplied are from multiple chromosomes")
 if (thechrs != goodchr) stop(paste("probeids requested not all on chr", chr))
 tss = mget( probeids, get(paste(gsub(".db", "", anno), "CHRLOC", sep="")), ifnotfound=NA)
 ends = mget( probeids, get(paste(gsub(".db", "", anno), "CHRLOCEND", sep="")), ifnotfound=NA)
 tss = sapply(tss, "[", 1)
 ends = sapply(ends, "[", 1)
 if (any(isn <- (is.na(tss) | is.na(ends)))) {
    tss = tss[-which(isn)]
    ends = ends[-which(isn)]
    probeids = probeids[-which(isn)]
 }
 if (!useEnd) ends = tss
 ans = GRanges(IRanges(start=abs(tss)-radius, end=abs(ends)+radius), seqnames=chr)
 names(ans) = probeids
 ans
}

snpIdsCisToGenes = function( mgr, chr, snpGR, radius=5e5, useEnd=FALSE ) {
 allgenes = colnames(mgr@fflist[[1]])
 CR = cisRanges(allgenes, chr=chr, anno=mgr@geneanno, radius=radius, useEnd=useEnd)
 FF = findOverlaps(snpGR, CR )
 allrs = names(snpGR)
 cisinds = split(FF@matchMatrix[,1], FF@matchMatrix[,2])
 cisrs = lapply(cisinds, function(x) allrs[x] )
 names(cisrs) = names(CR)[ FF@matchMatrix[,2][-which(duplicated(FF@matchMatrix[,2])) ]]
 cisrs
}

#OLDcisScores = function( mgr, ffind=1, chr, snpGR, radius=5e5, applier=lapply ) {
# cisrs = snpIdsCisToGenes( mgr, chr, snpGR, radius )
# onboard = rownames(mgr@fflist[[ffind]])
# cisrs = lapply(cisrs, function(x) intersect(onboard, x))
# ans = applier(1:length(cisrs), function(x) {
#      scores = as.ram(mgr@fflist[[ffind]][ cisrs[[x]], names(cisrs)[x] ] )/mgr@shortfac
#      names(scores) = cisrs[[x]]
#      scores
#      })
# names(ans) = names(cisrs)
# ans
#}

cisScores = function (mgr, ffind = 1, chr, snpGR, radius = 5e+05, applier = lapply, 
    minMAF = 0, minGTF = 0, useEnd=FALSE) 
{
    cisrs = snpIdsCisToGenes(mgr, chr, snpGR, radius, useEnd=useEnd)
    onboard = rownames(mgr@fflist[[ffind]])
    cisrs = lapply(cisrs, function(x) intersect(onboard, x))
    ans = applier(1:length(cisrs), function(x) {
        scores = as.ram(mgr@fflist[[ffind]][cisrs[[x]], names(cisrs)[x]])/mgr@shortfac
        names(scores) = cisrs[[x]]
        scores
    })
    okinds = NULL
    if (minMAF > 0) {
        maf = as.numeric(mgr@summaryList[[ffind]][, "MAF"])/mgr@shortfac
        okinds = which(maf >= minMAF)
    }
    if (minGTF > 0) {
        mgtf = as.numeric(mgr@summaryList[[ffind]][, "mGTF"])/mgr@shortfac
        tmp = which(mgtf >= minGTF)
        if (!is.null(okinds)) 
            okinds = intersect(tmp, okinds)
        else okinds = tmp
    }
    if (!is.null(okinds)) {
        okrs = rownames(mgr@summaryList[[ffind]])[okinds]
        ans = lapply(ans, function(x) x[intersect(okrs, names(x))])
    }
    names(ans) = names(cisrs)
    ans
}

eqtlTestsMACH = function(smlSet, machmat, rhs=~1-1,
   runname="foo", targdir="foo", geneApply=lapply, chromApply=lapply,
   shortfac = 100, computeZ=FALSE, family, ... ) {
 theCall = match.call()
 sess = sessionInfo()
 if (missing(family)) family="gaussian"
 fnhead = paste(targdir, "/", runname, "_", sep="")
 geneNames = featureNames(smlSet)
 chrNames = names(smList(smlSet))
 ngenes = length(geneNames)
 nchr = length(chrNames)
 system(paste("mkdir", targdir))
 cres = chromApply( chrNames, function(chr) {
   snpdata = machmat # smList(smlSet)[[chr]]
   #targff = paste( fnhead, "chr", chr, "_", "g", gene, ".ff" , sep="" )
   targff = paste( fnhead, "chr", chr, ".ff" , sep="" )
   snpnames = rownames(snpdata) #colnames(snpdata)
   nsnps = nrow(snpdata) #ncol(snpdata)
   store = ff( initdata=0, dim=c(nsnps, ngenes), dimnames=list(snpnames, geneNames), vmode="short",
                 filename = targff )
   geneApply( geneNames, function(gene) {
     ex = exprs(smlSet)[gene,]
     fmla = formula(paste("ex", paste(as.character(rhs),collapse=""), collapse=" "))
     numans = snp.rhs.testsMACH(fmla, snp.data=snpdata, data=pData(smlSet), family=family, ...)@chisq
     if (computeZ) {
       stop("not handled")
       numans = sqrt(numans)
       signl = snp.rhs.estimates( fmla, snp.data=snpdata, data=pData(smlSet), family="gaussian", ... )
       bad = which(unlist(lapply(signl, is.null)))
       if (length(bad)>0) signl[bad] = list(beta=NA)
       ifelse(unlist(signl)>=0, 1, -1)
       numans = numans*signl
     }
     miss = is.na(numans)
     if (any(miss) & !computeZ) numans[which(miss)] = rchisq(length(which(miss)), 1)
     if (any(miss) & computeZ) numans[which(miss)] = rnorm(length(which(miss)))
     store[, gene, add=TRUE] = shortfac*numans
     NULL
     }) # end gene apply
  store
  })  # end chr apply
  names(cres) = chrNames
  exdate = date()
  new("eqtlTestsManager", fflist=cres, call=theCall, sess=sess, 
        exdate=exdate, shortfac=shortfac, geneanno=annotation(smlSet),
        df=1)
}

manhPlot = function( probeid, mgr, ffind, namedlocvec=NULL, locGRanges=NULL,
   plotter=smoothScatter, tx=function(x)-log10(1-pchisq(x,1)), 
   xlab = paste("pos. on ",names(fflist(mgr))[ffind]),
   ylab = "-log10 p", ... ) {
 if (!(is(mgr, "eqtlTestsManager"))) stop("mgr must inherit from eqtlTestsManager")
 if (is.null(namedlocvec) & is.null(locGRanges)) stop("one of namedlocvec and locGRanges must be non-null")
 if (is.null(namedlocvec) & is.null(names(locGRanges))) stop("locGRanges must have non-null names")
 vals = mgr[,  probeId(probeid), drop=FALSE]
 vals = vals[[ffind]][,]
 rsidInVals = names(vals)
 if (!is.null(locGRanges)) {
   rsidInLocs = names(locGRanges)
   namedlocvec = start(locGRanges)
   names(namedlocvec) = names(locGRanges)
   }
 okrs = intersect(rsidInVals, names(namedlocvec))
 mm = match(okrs, names(namedlocvec))
 vv = match(okrs, names(vals))
 loc = namedlocvec[mm]
 vals = as.numeric(vals[vv])
 plotter(loc, tx(vals), xlab=xlab, ylab=ylab, ...)
 anno = mgr@geneanno
 if (require(anno, character.only=TRUE)) {
   packref = function(tag="CHRLOC") get(paste(gsub(".db", "", anno), tag, sep=""))
   gloc = get(probeid, packref())
   axis(3, label=get(probeid, packref("SYMBOL")),
           at=abs(gloc[1]), col="red", lwd=2)
   }
 invisible(NULL)
}
 
meqtlTests = function(listOfSmls, rhslist,
   runname="mfoo", targdir="mfoo", geneApply=lapply, chromApply=lapply,
   shortfac = 100, computeZ=FALSE, harmonizeSNPs = FALSE, uncert=TRUE, 
   saveSummaries=TRUE, family, ... ) {
 theCall = match.call()
 sess = sessionInfo()
 if (missing(family)) family="gaussian"
 allfeat = lapply(listOfSmls, featureNames)
 smlSet1 = listOfSmls[[1]]
 fint = allfeat[[1]]
 for (i in 2:length(allfeat)) fint = intersect(fint, allfeat[[i]])
 if (length(fint)==0) stop("null intersection of featureNames for smlSet list elements")
 listOfSmls = reduceGenes( listOfSmls, probeId(fint) )
 if (harmonizeSNPs) listOfSmls = makeCommonSNPs( listOfSmls )
  else if(!isTRUE(checkCommonSNPs( listOfSmls ))) stop("harmonizeSNPs = FALSE but SNPs not common across listOfSmls, run makeCommonSNPs")

 smlSet1 = listOfSmls[[1]]
 fnhead = paste(targdir, "/", runname, "_", sep="")
 geneNames = featureNames(smlSet1)
 chrNames = names(smList(smlSet1))
 ngenes = length(geneNames)
 nchr = length(chrNames)
 system(paste("mkdir", targdir))
#
# there will be one ff file per chromosome which will accumulate
# all information across smlSets
#
 targffs = paste( fnhead, "chr", chrNames, ".ff", sep="" )
 allSnpnames = lapply(smList(listOfSmls[[1]]), colnames)
 ffRefList = lapply( 1:nchr, function(chr)
    ff( initdata = 0, dim=c( length(allSnpnames[[chr]]), ngenes),
        dimnames = list(allSnpnames[[chr]], geneNames), vmode="short",
        filename=targffs[chr] ))
 names(ffRefList) = chrNames
 
# chopped from eqtlTests
 summfflist = list()
 if (saveSummaries) {
  # get MAF and minGTF for all SNP
  sumfn = paste(fnhead, chrNames, "_summ.ff", sep="")
  if ("multicore" %in% search()) {
    summfflist = mclapply( 1:length(chrNames), function(i) ffSnpSummary(smList(smlSet)[[i]], sumfn[i],
         fac=shortfac))
    } else {
          for (i in 1:length(sumfn))
              summfflist[[chrNames[i]]] = ffSnpSummary(smList(smlSet)[[i]], sumfn[i])
          }
  # ok, now just save references in object
  }

 cres = chromApply( chrNames, function(chr) {
  for (theSS in 1:length(listOfSmls)) {
   smlSet = listOfSmls[[theSS]]
   store = ffRefList[[chr]]
   snpdata = smList(smlSet)[[chr]]
   geneApply( geneNames, function(gene) {
     ex = exprs(smlSet)[gene,]
     fmla = formula(paste("ex", paste(as.character(rhslist[[theSS]]),collapse=""), collapse=" "))
     numans = snp.rhs.tests(fmla, snp.data=snpdata, 
         data=pData(smlSet), family=family, uncertain=uncert, ...)@chisq
     if (computeZ) {
       numans = sqrt(numans)
       signl = snp.rhs.estimates( fmla, snp.data=snpdata, data=pData(smlSet), family=family, ... )
       bad = which(unlist(lapply(signl, is.null)))
       if (length(bad)>0) signl[bad] = list(beta=NA)
       ifelse(unlist(signl)>=0, 1, -1)
       numans = numans*signl
     }
     miss = is.na(numans)
     if (any(miss) & !computeZ) numans[which(miss)] = rchisq(length(which(miss)), 1)
     if (any(miss) & computeZ) numans[which(miss)] = rnorm(length(which(miss)))
     store[, gene, add=TRUE] = shortfac*numans
     NULL
     }) # end gene apply
   } # end iterate over smlSet list
   store
  })  # end chr apply
  names(cres) = chrNames
  exdate = date()
  new("eqtlTestsManager", fflist=cres, call=theCall, sess=sess, 
        exdate=exdate, shortfac=shortfac, geneanno=annotation(smlSet1),
        df=length(listOfSmls), summaryList=summfflist)
}
