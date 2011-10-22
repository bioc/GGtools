
restrictProbesToChrom = function(smlSet, chrom) {
 ganno = smlSet@annotation
 require(ganno, character.only=TRUE)
 annomappref = gsub(".db", "", ganno)
 map = get(mapn <- paste(annomappref, "CHR", sep=""))
 if (!chrom %in% mappedkeys(revmap(map))) stop(paste(chrom, "not mapped in", mapn))
 psn = get(chrom, revmap(map))
 ans = smlSet[ probeId(intersect(psn, featureNames(smlSet))), ]
 if (nrow(ans) == 0) stop("no probes present after filter")
 ans
}

.restrictProbesToChromGRL = function(smlSet, chrom, GRL) {
#
# at this point it is assumed that chrom is a single integer or string representing an integer
#
 chrtok = paste("chr", chrom, sep="")
 cn = as(seqnames(GRL), "character")
 if (!(chrtok %in% cn)) stop(paste("can't find", chrtok, "in seqnames(GRL)"))
 probesOnChr = names(GRL)[ which(cn == chrtok) ]
 if (length(probesOnChr) < 1) stop(paste("no probes on chr", chrtok, "found in geneGRL"))
 okpr = intersect(probesOnChr, featureNames(smlSet))
 if (length(okpr) < 1) stop("no probes on chr", chrtok, "resident in both geneGRL and featureNames(smlSet)")
 smlSet[ probeId(okpr), ]
}

restrictProbesToChromGRL = function (smlSet, chrom, GRL)
{
    chrtok = paste("chr", chrom, sep = "")
    GR = GRL[[chrtok]]
    cn = as(seqnames(GR), "character")
    if (!(chrtok %in% cn))
        stop(paste("can't find", chrtok, "in seqnames(GRL)[[chrom]]"))
    probesOnChr = names(GR)[which(cn == chrtok)]
    if (length(probesOnChr) < 1)
        stop(paste("no probes on chr", chrtok, "found in geneGRL"))
    okpr = intersect(probesOnChr, featureNames(smlSet))
    if (length(okpr) < 1)
        stop("no probes on chr", chrtok, "resident in both geneGRL and featureNames(smlSet)")
    smlSet[probeId(okpr), ]
}


setClass("multiCisDirector", representation(
  mgrs="list"))
setMethod("show", "multiCisDirector", function(object){
 cat("multiCisDirector instance with", length(object@mgrs), "eqtlTestsManagers.\n")
 cat("elements:", selectSome(names(object@mgrs)), "\n")
})

getbds = function (mat) 
{
#
# produces char vector of form 1K.5K from matrix rows 1000, 5000, for example
#
    bds = cbind(mat, apply(mat, 1, sum))[, c(1, 3),drop=FALSE] # key to set drop for single radius
    fm = gsub(" ", "", format(bds))
    fm = gsub("000000$", "M", fm)
    fm = gsub("000$", "K", fm)
    apply(fm, 1, paste, collapse = ".")
}


collectGeneRanges = function(mcd, chrpref="chr", applier=lapply) {
 if (!is(mcd, "multiCisDirector")) stop("needs multiCisDirector instance")
 chrs = names(mcd@mgrs)
 if (is.null(chrs)) stop("managers in mcd must be named")
 nmgrs = length(mcd@mgrs)
 ganno = mcd@mgrs[[1]]@geneanno   # assume
 require(ganno, character.only=TRUE)
 ganno = gsub(".db", "", ganno)
 gloc = get(paste(ganno, "CHRLOC", sep=""))
 glocend = get(paste(ganno, "CHRLOCEND", sep=""))
 allplist = lapply( mcd@mgrs, probeNames )
 stlist = lapply( allplist, function(x) abs(sapply(AnnotationDbi::mget(x, gloc), "[", 1)))
 enlist = lapply( allplist, function(x) abs(sapply(AnnotationDbi::mget(x, glocend), "[", 1)))
 oklist = lapply( 1:nmgrs, function(x) 
              which(!is.na(stlist[[x]]) & !is.na(enlist[[x]])))
 stlist = lapply( 1:nmgrs, function(i) stlist[[i]][ oklist[[i]] ] )
 enlist = lapply( 1:nmgrs, function(i) enlist[[i]][ oklist[[i]] ] )
 allplist = lapply( 1:nmgrs, function(i) allplist[[i]][ oklist[[i]] ] )
 require(GenomicRanges)
 ans = lapply( 1:nmgrs, function(x) {
     tmp = GRanges(seqnames=chrs[x], IRanges(stlist[[x]], enlist[[x]]) ) 
     names(tmp) = allplist[[x]]
     tmp
     } )
 names(ans) = chrs
 ans
}


collectSNPRanges = function(mcd, 
     slpref="ch", sprefInMgr="chr", applier=lapply,
     snpannopack) {
 if (!is(mcd, "multiCisDirector")) stop("needs multiCisDirector instance")
 chrs = names(mcd@mgrs)
 if (is.null(chrs)) stop("managers in mcd must be named")
 nmgrs = length(mcd@mgrs)
 chrForSNPlocs = gsub(sprefInMgr, slpref, chrs)
 if (length(grep(slpref, chrForSNPlocs))==0) # guess only numeric chrnames used
   chrForSNPlocs = paste(slpref, chrForSNPlocs, sep="")
 require(snpannopack, character.only=TRUE)
 allsr = applier( 1:nmgrs, function(i) {
  tmpsl = getSNPlocs(chrForSNPlocs[i])
  rsnames = paste("rs", tmpsl$RefSNP_id, sep="")
  rsinmgr = rownames(mcd@mgrs[[i]]@fflist[[1]])
  tmpr = GRanges(seqnames=chrs[i], IRanges(tmpsl$loc, width=1))
  names(tmpr) = rsnames
  tmpr[ intersect( rsnames, rsinmgr ) ]
  })
 names(allsr) = chrs
 allsr
}

setClass("cisProxScores", representation(call="call"), contains="list")
setMethod("show", "cisProxScores", function(object) {
cat("GGtools cisProxScores instance.\n")
cat("The call was: ")
print(object@call)
cat("intervals examined:", selectSome(names(object)), "\n")
})


cisProxScores = function( smlSet, fmla, dradset, direc=NULL,
   folder, runname, geneApply=lapply, saveDirector=TRUE, 
   snpGRL = NULL, geneGRL = NULL, 
   snpannopack="SNPlocs.Hsapiens.dbSNP.20100427", ffind=NULL,
   ... ) {
  thecall = match.call()
  if (is.null(ffind)) stop("must set ffind (usually to 1)")
  if (is.null(direc)) {
   chrs = names(smList(smlSet))
   nchrs = gsub("chr", "", chrs)
   mgrs = lapply( 1:length(chrs), function(c) {
    if (is.null(geneGRL)) {
     thisset = restrictProbesToChrom( smlSet, nchrs[c] )
     } else {
     thisset = restrictProbesToChromGRL( smlSet, nchrs[c], geneGRL)
     }
     thisset = thisset[ chrnum(chrs[c]), ]
     eqtlTests( thisset, fmla, targdir=folder, runname=paste(runname, "_",
        c, sep=""), geneApply=geneApply, ... )
     } )
   names(mgrs) = chrs
   direc = new("multiCisDirector", mgrs=mgrs)
   if (saveDirector) {
        dirn = paste(folder, "_director", sep="")
        assign(dirn, direc)
        save(list=dirn, file=paste(dirn, ".rda", sep=""))
        }
  }
  if (is.null(geneGRL)) {
   gr = collectGeneRanges(direc, applier=geneApply) # reuse of geneApply not ideal
    } else gr = geneGRL
  if (!isTRUE(all(names(direc@mgrs) %in% names(gr))))
         stop("geneGRL must be a list of GRanges with all names including names(direc@mgrs)")
  if (is.null(snpGRL)) {
   sr = collectSNPRanges(direc, applier=geneApply, snpannopack=snpannopack)
   } else sr = snpGRL
  if (!isTRUE(all(names(direc@mgrs) %in% names(sr))))
         stop("snpGRL must be a list of GRanges with names including names(direc@mgrs)")
  
  if (any(diff(dradset)<0)) stop("diff(dradset) must yield only positive numbers")
  radmat = cbind(c(0, dradset[-length(dradset)]), c(dradset[1], diff(dradset)))
  radnms = paste("FL", getbds(radmat), sep="")
 # following code will put gene+dradset hole in first element... not intended, see fill below
  intlist = lapply(1:nrow(radmat), function(z) {  # over family of radii
    ans = lapply(gr, function(gra) {  # over chrom-specific granges
        curg = gra
        ranges(curg) = ranges(curg)+radmat[z,1]
        flankingOnly(curg, radmat[z,2]) 
    } ) 
   ans 
  } )
  intlist[[1]] = lapply(gr, function(gra) {
       curg = gra
       ranges(curg) = ranges(curg)+dradset[1]
       curg
       }
      )  # fill initial hole
  names(intlist) = radnms
  if (TRUE) {        #  how can you reuse "gr" below?  it seems right because gr has been
                     #  transformed above to intlist but potentially confusing
  ## there is a problem with the original coding if the order of spaces
  ## in managers and address interval lists is not common.  for example the
  ## intervals might have names chr1 chr10 chr11 and so on...
  ## need to iterate over the names and select directly
  allspaces = intersect( names(direc@mgrs), names(gr) )
  allspaces = intersect( allspaces, names(sr) )
  ans = lapply( intlist, function(ingr) {   # recycled "gr" here previously
     cans = lapply( allspaces, function(mgrind) {
      cat(mgrind)
      scoresInRanges( direc@mgrs[[mgrind]], ingr[[mgrind]], sr[[mgrind]],
        applier=geneApply, ffind=ffind ) } ) 
     names(cans) = allspaces # names(direc@mgrs)
     cans
     } 
    )
  } 
  new("cisProxScores", ans, call=thecall)  # final value
}


#mcisProxScores = function( listOfSmlSets, listOfFmlas, dradset, direc=NULL,
#   folder, runname, geneApply=mclapply, saveDirector=TRUE, 
#   makeCommonSNPs=FALSE,
#   snpGRL = NULL, geneGRL = NULL,
#   snpannopack="SNPlocs.Hsapiens.dbSNP.20100427", ffind=NULL,
#   ... ) {
#  if (is.null(ffind)) stop("must set ffind (usually to 1)")
#  if (missing(dradset)) stop("must supply dradset")
#  thecall = match.call()
#  if (is.null(direc)) {
#  if (length(listOfSmlSets) < 2) stop("need list of > 1 smlSet")
#  if (makeCommonSNPs) listOfSmlSets = makeCommonSNPs(listOfSmlSets)
#  fnl = lapply(listOfSmlSets, featureNames)
#  fn1 = fnl[[1]]
#  for (i in 2:length(fnl)) if (!all.equal(fnl[[i]], fn1)) stop("need congruent featureSets for all input smlSets")
#    }
#  if (is.null(direc)) {
#   chrs = names(smList(listOfSmlSets[[1]]))
#   nchrs = gsub("chr", "", chrs)
#   mgrs = lapply( 1:length(chrs), function(c) {
#    if (is.null(geneGRL)) {
#     thislist = lapply(listOfSmlSets, function(s) restrictProbesToChrom(
#         s, nchrs[c]))
#     } else {
#     thislist = lapply(listOfSmlSets, function(s) restrictProbesToChromGRL(
#         s, nchrs[c], geneGRL))
#     }
#     thislist = lapply(thislist, function(x) x[chrnum(chrs[c]),] )
#     meqtlTests( thislist, listOfFmlas, targdir=folder, runname=paste(runname, "_",
#        c, sep=""), geneApply=geneApply, ... )
#     } )
#   names(mgrs) = chrs
#   direc = new("multiCisDirector", mgrs=mgrs)
#   if (saveDirector) {
#        dirn = paste(folder, "_director", sep="")
#        assign(dirn, direc)
#        save(list=dirn, file=paste(dirn, ".rda", sep=""))
#        }
#  }
#  if (is.null(geneGRL)) {
#   gr = collectGeneRanges(direc, applier=geneApply) # reuse of geneApply not ideal
#    } else gr = geneGRL
#  if (!isTRUE(all(names(direc@mgrs) %in% names(gr))))
#         stop("geneGRL must be a list of GRanges with all names including names(direc@mgrs)")
#  if (is.null(snpGRL)) {
#   sr = collectSNPRanges(direc, applier=geneApply, snpannopack=snpannopack)
#   } else sr = snpGRL
#  if (!isTRUE( all(names(direc@mgrs) %in% names(sr))))
#         stop("snpGRL must be a list of GRanges with names including names(direc@mgrs)")
#  
#  if (any(diff(dradset)<0)) stop("diff(dradset) must yield only positive numbers")
#  radmat = cbind(c(0, dradset[-length(dradset)]), c(dradset[1], diff(dradset)))
#  radnms = paste("FL", getbds(radmat), sep="")
#  intlist = lapply(1:nrow(radmat), function(z) {  # over family of radii
#    ans = lapply(gr, function(gra)  # over chrom-specific granges
#        flankingOnly(gra+radmat[z,1], radmat[z,2]) ) 
#    ans
#    } )
#  intlist[[1]] = lapply(gr, function(gra) gra+dradset[1])  # fill initial hole
#  names(intlist) = radnms
#  allspaces = intersect( names(direc@mgrs), names(gr) )
#  allspaces = intersect( allspaces, names(sr) )
#  if (TRUE) {
#  ans = lapply( intlist, function(gr) {  # this binding of gr is annoying...
#     cans = lapply( allspaces, function(mgrind) {
#      cat(mgrind)
#      scoresInRanges( direc@mgrs[[mgrind]], gr[[mgrind]], sr[[mgrind]],
#        applier=geneApply, ffind=ffind ) } ) 
#     names(cans) = allspaces # names(direc@mgrs)
#     cans
#     } 
#    )
#  } 
#  new("cisProxScores", ans, call=thecall)  # final value
#}
#


#fmcisRun = function( listOfSmlSets, listOfFmlas, 
#   folder, runname, geneApply=mclapply, saveDirector=TRUE, ffind=1,
#   geneGRL, snpGRL )
#    {
#   chrs = names(smList(listOfSmlSets[[1]]))
#   nchrs = gsub("chr", "", chrs)
#   mgrs = lapply( 1:length(chrs), function(c) {
#     thislist = lapply(listOfSmlSets, function(s) restrictProbesToChromGRL(
#         s, nchrs[c], geneGRL))
#     thislist = lapply(thislist, function(x) x[chrnum(chrs[c]),] )
#     mfeqtltests( listOfEgtSet=thislist, listOfRhs=listOfFmlas, nchunk=10, targdir=folder, runname=paste(runname, "_",
#        c, sep=""), geneApply=geneApply )
#     } )
#   names(mgrs) = chrs
#   direc = new("multiCisDirector", mgrs=mgrs)
#   if (saveDirector) {
#        dirn = paste(folder, "_director", sep="")
#        assign(dirn, direc)
#        save(list=dirn, file=paste(dirn, ".rda", sep=""))
#        }
#   direc
#  }



scoresByGenes = function (cps, intvind = 1, as.GRanges=TRUE, dups2max=TRUE, snpGR=NULL,
  scoreConverter=function(x)x ) { 
#
# method for working with a cisProxScores result, to retrieve scores conveniently
#
    sco = lapply(cps[[intvind]], lapply, function(x) {
        tmp = x[[1]][, 1]
        names(tmp) = rownames(x[[1]])
        tmp
    })[[1]]
    if (!as.GRanges) return(lapply(sco,scoreConverter))
    if (is.null(snpGR)) stop("for as.GRanges, please supply snpGR with GRanges instance of SNP addresses")
    if (is.null(names(snpGR))) stop("snpGR names must hold the dbSNP id or other snp identifier")
    if (!dups2max) stop("no duplicate policy other than dups2max at present... please set dups2max to TRUE")
    allsco = unlist(sco)
    alln = unlist(lapply(sco, names))
    getd = split(allsco, alln)
    maxs = sapply(getd,max)
    names(maxs) = names(getd)
    okn = intersect(names(maxs), names(snpGR))
    maxs = maxs[okn]
    okn = match(names(maxs), names(snpGR), nomatch=0)
    snpGR = snpGR[okn[okn>0]]
    if (!all.equal(names(snpGR), names(maxs))) stop("something is wrong in match of snpGR names to max scores...")
    elementMetadata(snpGR)$score = scoreConverter(maxs)
    snpGR
}

