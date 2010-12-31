
restrictProbesToChrom = function(smlSet, chrom) {
 ganno = smlSet@annotation
 require(ganno, character.only=TRUE)
 annomappref = gsub(".db", "", ganno)
 map = get(mapn <- paste(annomappref, "CHR", sep=""))
 if (!chrom %in% mappedkeys(revmap(map))) stop(paste(chr, "not mapped in", mapn))
 psn = get(chrom, revmap(map))
 ans = smlSet[ probeId(intersect(psn, featureNames(smlSet))), ]
 if (nrow(ans) == 0) stop("no probes present after filter")
 ans
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
 stlist = lapply( allplist, function(x) abs(sapply(mget(x, gloc), "[", 1)))
 enlist = lapply( allplist, function(x) abs(sapply(mget(x, glocend), "[", 1)))
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


collectSNPRanges = function(mcd, slpref="ch", sprefInMgr="chr", applier=lapply) {
 if (!is(mcd, "multiCisDirector")) stop("needs multiCisDirector instance")
 chrs = names(mcd@mgrs)
 if (is.null(chrs)) stop("managers in mcd must be named")
 nmgrs = length(mcd@mgrs)
 chrForSNPlocs = gsub(sprefInMgr, slpref, chrs)
 if (length(grep(slpref, chrForSNPlocs))==0) # guess only numeric chrnames used
   chrForSNPlocs = paste(slpref, chrForSNPlocs, sep="")
 require(SNPlocs.Hsapiens.dbSNP.20100427)
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


.cisProxScores = function( smlSet, fmla, dradset, direc=NULL,
   folder, runname, geneApply=mclapply, saveDirector=TRUE, 
   geneCentric = TRUE, retain=10, ... ) {
  thecall = match.call()
  if (is.null(direc)) {
   chrs = names(smList(smlSet))
   nchrs = gsub("chr", "", chrs)
   mgrs = lapply( 1:length(chrs), function(c) {
     thisset = restrictProbesToChrom( smlSet, nchrs[c] )
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
  gr = collectGeneRanges(direc, applier=geneApply) # reuse of geneApply not ideal
  sr = collectSNPRanges(direc, applier=geneApply)
  if (any(diff(dradset)<0)) stop("diff(dradset) must yield only positive numbers")
  radmat = cbind(c(0, dradset[-length(dradset)]), c(dradset[1], diff(dradset)))
  radnms = paste("FL", getbds(radmat), sep="")
  intlist = lapply(1:nrow(radmat), function(z) {  # over family of radii
    ans = lapply(gr, function(gra)  # over chrom-specific granges
        flankingOnly(gra+radmat[z,1], radmat[z,2]) ) 
    ans
    } )
  names(intlist) = radnms
  if (geneCentric) {
  ans = lapply( intlist, function(gr) {
     cans = lapply( 1:length(direc@mgrs), function(mgrind) {
      cat(mgrind)
      scoresInRanges( direc@mgrs[[mgrind]], gr[[mgrind]], sr[[mgrind]],
        applier=geneApply ) } ) 
     names(cans) = names(direc@mgrs)
     cans
     } 
    )
  } else {  # end geneCentric
  # for snpcentric reporting, filter SNP to proximity ranges
    srtargs2 = lapply(1:length(intlist), function(fammem) {
         ans = lapply(1:length(sr), function(chr) sr[[chr]][
             which(IRanges::"%in%"(ranges(sr[[chr]]), ranges(intlist[[fammem]][[chr]]))) ] )
         names(ans) = names(sr)
         ans
         } )
    names(srtargs2) = names(intlist)
  # now collect the topFeats results for each SNP for a specified number of genes
  #  per SNP
    ans = lapply(1:length(srtargs2), function(fammem) {
      ans = lapply(1:length(srtargs2[[fammem]]), function(chr) {
        rs = names(srtargs2[[fammem]][[chr]])
        ans = lapply(rs, function(sn) {
              topFeats(rsid(sn), mgr=direc@mgrs[[chr]],
         ffind=1, useSym=FALSE, n=retain ) }  )
        names(ans) = rs
        ans
       } ) # done chr
      names(ans) = names(srtargs2[[fammem]])
      ans
      } ) # done fammem
    names(ans) = names(srtargs2)
    }  # conclude snp-centric chunk
   new("cisProxScores", ans, call=thecall)  # final value
}
    
cisProxScores = function( smlSet, fmla, dradset, direc=NULL,
   folder, runname, geneApply=mclapply, saveDirector=TRUE, 
   ... ) {
  thecall = match.call()
  if (is.null(direc)) {
   chrs = names(smList(smlSet))
   nchrs = gsub("chr", "", chrs)
   mgrs = lapply( 1:length(chrs), function(c) {
     thisset = restrictProbesToChrom( smlSet, nchrs[c] )
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
  gr = collectGeneRanges(direc, applier=geneApply) # reuse of geneApply not ideal
  sr = collectSNPRanges(direc, applier=geneApply)
  if (any(diff(dradset)<0)) stop("diff(dradset) must yield only positive numbers")
  radmat = cbind(c(0, dradset[-length(dradset)]), c(dradset[1], diff(dradset)))
  radnms = paste("FL", getbds(radmat), sep="")
  intlist = lapply(1:nrow(radmat), function(z) {  # over family of radii
    ans = lapply(gr, function(gra)  # over chrom-specific granges
        flankingOnly(gra+radmat[z,1], radmat[z,2]) ) 
    ans
    } )
  names(intlist) = radnms
  if (TRUE) {
  ans = lapply( intlist, function(gr) {
     cans = lapply( 1:length(direc@mgrs), function(mgrind) {
      cat(mgrind)
      scoresInRanges( direc@mgrs[[mgrind]], gr[[mgrind]], sr[[mgrind]],
        applier=geneApply ) } ) 
     names(cans) = names(direc@mgrs)
     cans
     } 
    )
  } 
  new("cisProxScores", ans, call=thecall)  # final value
}


mcisProxScores = function( listOfSmlSets, listOfFmlas, dradset, direc=NULL,
   folder, runname, geneApply=mclapply, saveDirector=TRUE, 
   makeCommonSNPs=FALSE,
   ... ) {
  thecall = match.call()
  if (length(listOfSmlSets) < 2) stop("need list of > 1 smlSet")
  if (makeCommonSNPs) listOfSmlSets = makeCommonSNPs(listOfSmlSets)
  fnl = lapply(listOfSmlSets, featureNames)
  fn1 = fnl[[1]]
  for (i in 2:length(fnl)) if (!all.equal(fnl[[i]], fn1)) stop("need congruent featureSets for all input smlSets")
  if (is.null(direc)) {
   chrs = names(smList(listOfSmlSets[[1]]))
   nchrs = gsub("chr", "", chrs)
   mgrs = lapply( 1:length(chrs), function(c) {
     thislist = lapply(listOfSmlSets, function(s) restrictProbesToChrom(
         s, nchrs[c]))
     thislist = lapply(thislist, function(x) x[chrnum(chrs[c]),] )
     meqtlTests( thislist, listOfFmlas, targdir=folder, runname=paste(runname, "_",
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
  gr = collectGeneRanges(direc, applier=geneApply) # reuse of geneApply not ideal
  sr = collectSNPRanges(direc, applier=geneApply)
  if (any(diff(dradset)<0)) stop("diff(dradset) must yield only positive numbers")
  radmat = cbind(c(0, dradset[-length(dradset)]), c(dradset[1], diff(dradset)))
  radnms = paste("FL", getbds(radmat), sep="")
  intlist = lapply(1:nrow(radmat), function(z) {  # over family of radii
    ans = lapply(gr, function(gra)  # over chrom-specific granges
        flankingOnly(gra+radmat[z,1], radmat[z,2]) ) 
    ans
    } )
  names(intlist) = radnms
  if (TRUE) {
  ans = lapply( intlist, function(gr) {
     cans = lapply( 1:length(direc@mgrs), function(mgrind) {
      cat(mgrind)
      scoresInRanges( direc@mgrs[[mgrind]], gr[[mgrind]], sr[[mgrind]],
        applier=geneApply ) } ) 
     names(cans) = names(direc@mgrs)
     cans
     } 
    )
  } 
  new("cisProxScores", ans, call=thecall)  # final value
}
