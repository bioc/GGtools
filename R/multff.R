
makeCommonSNPs = function( listOfSms ) {
  rsidlist = intersectSnps( listOfSms )
  trimSnps( listOfSms, rsidlist )
}

rsnum = function(x) lapply(x$fflist, rownames)

 
multffCT = function(listOfSms, gfmlaList, geneinds=1:10, harmonizeSNPs=FALSE, 
     targdir=".", runname="foo", overwriteFF=TRUE, fillNA=TRUE, ncores=2, mc.set.seed=TRUE, vmode="single", shortfac=100, ...) {
  theCall = match.call()
  if (!file.exists(targdir)) stop("targdir must exist prior to invocation of multffCT")
  require(ff, quietly=TRUE)
  .checkArgsMF( listOfSms, gfmlaList, geneinds, targdir, runname )
  listOfSms = reduceGenes( listOfSms, geneinds )
  if (harmonizeSNPs) listOfSms = makeCommonSNPs( listOfSms )
  else if(!isTRUE(checkCommonSNPs( listOfSms ))) stop("harmonizeSNPs = FALSE but SNPs not common across listOfSms, run makeCommonSNPs")
  sumScores2ff( listOfSms, gfmlaList, targdir, runname, theCall, overwriteFF=overwriteFF,
       fillNA=fillNA, write=TRUE, ncores=ncores, vmode=vmode, mc.set.seed=mc.set.seed, shortfac=shortfac, ... )
}


.checkArgsMF = function( listOfSms, gfmlaList, geneinds, targdir, runname ) {
  if (!inherits(listOfSms, "list")) stop("listOfSms must inherit from list")
  if (length(listOfSms) != length(gfmlaList)) stop("gfmlaList must have same length as listOfSms")
  allc = sapply(listOfSms, function(x) inherits(x, "smlSet"))
  if (!isTRUE(all(allc))) stop("each element of listOfSms must inherit from GGtools smlSet")
 allcsets = sapply(listOfSms, function(x) names(smList(x)))
 if(!is.atomic(allcsets)) stop("probably not all elements of listOfSms have same chromosome set")
 allfn = unique(unlist(fnlist <- lapply(listOfSms, featureNames)))
 if (length(fnlist)>1) {
  for (i in 2:length(fnlist)) 
     if (!isTRUE(all.equal(fnlist[[i]], fnlist[[1]]))) stop("all smlSets must have same feature set")
 }
 allfi = sapply(lapply(listOfSms, featureNames), length)
 if (inherits(geneinds, "character") & !(all(geneinds %in% allfn)))
	stop("some geneinds not in featureNames of listOfSms elements")
 else if (inherits(geneinds, "numeric") & !(max(geneinds)<=max(allfi)))
	stop("some geneinds not in 1:length(featureNames(sms)) for some sms in listOfSms")
 }

reduceGenes = function( listOfSms, geneinds )
  lapply( listOfSms, function(x) x[ geneinds, ] )

intersectSnps = function( listOfSms ) {
  nsms = length(listOfSms)
  nchr = length(smList(listOfSms[[1]]))
  chnames = names(smList(listOfSms[[1]]))
  rsidlist = lapply(smList(listOfSms[[1]]), colnames)
  if (length(listOfSms) == 1) return(rsidlist)
  for (j in 2:nsms)
    for (k in 1:length(rsidlist))
      rsidlist[[k]]  = intersect(rsidlist[[k]], colnames(smList(listOfSms[[j]])[[k]]))
  rsidlist
}

trimSnps = function( listOfSms, rsidlist ) {
  ans = list()
  cat("making smlSets with common SNPs: \n")
  for (j in 1:length(listOfSms)) {
    cat("smlSet", j, "\nchrom ")
    tmp = listOfSms[[j]]
    for (k in 1:length(rsidlist)) {
      cat(k)
      tmp@smlEnv$smList[[k]] = listOfSms[[j]]@smlEnv$smList[[k]][, rsidlist[[k]]]
      }
    ans[[j]] = tmp
    cat("\n")
    }
  names(ans) = names(listOfSms)
  cat("done.\n")
  gc()
  ans
}

checkCommonSNPs = function( listOfSms ) {
  if (length(listOfSms) == 1) return(TRUE)
  rsidlist = lapply( smList(listOfSms[[1]]), colnames )
  for (j in 2:length(listOfSms) ) {
    tmp = lapply(smList(listOfSms[[j]]), colnames )
    if (length(tmp) != length(rsidlist)) stop("num chroms not common to all elements of listOfSms")
    for (k in 1:length(rsidlist))
      if (!isTRUE(all.equal(tmp[[k]], rsidlist[[k]]))) {
              cat("sms ", j, " chr ", k, "\n")
              stop("rsid lists not common across smlsets")
              }
     }
  return(TRUE)
}

sumScores2ff = function( listOfSms, gfmlaList, targdir, runname, theCall=call("1"), 
      overwriteFF=FALSE, fillNA=TRUE, write=TRUE, ncores, vmode, mc.set.seed, shortfac, ... ) {
  fnhead = paste(targdir, "/", runname, "_", sep="")
  expd = lapply(listOfSms, experimentData)
  nsms = length(listOfSms)
  nchr = length(smList(listOfSms[[1]]))
  nsnps = sapply(smList(listOfSms[[1]]), ncol)
  chrnames = names(smList(listOfSms[[1]])) 
  ngenes = length(featn <- featureNames(listOfSms[[1]]))
  generangetag = paste(featn[1], "_",  featn[length(featn)],sep="")
  gnames = featureNames(listOfSms[[1]])
  rslist = lapply(smList(listOfSms[[1]]), colnames)
  filenames = lapply( chrnames, function(x) paste(fnhead, "snpsOnChr", 
      x, "_", generangetag, "_df", nsms, ".ff", sep=""))
  if (write) fflist = lapply(1:nchr, function(x) ff(initdata=0.0, 
    dim=c(nsnps[x], ngenes), 
    dimnames = list(rslist[[x]], gnames),
    vmode=vmode,
    filename=filenames[[x]],
    overwrite=overwriteFF))
  else  fflist = lapply(1:nchr, function(x) ff(initdata=0.0, 
    dim=c(nsnps[x], ngenes), 
    dimnames = list(rslist[[x]], gnames),
    vmode=vmode))
#    filename=paste(fnhead, "snpsOnChr", x, "_", generangetag, ".ff", sep=""),
#    overwrite=overwriteFF)
#  )
  indmat = data.matrix(expand.grid( 1:nchr, 1:nsms ))[,c(2,1),drop=FALSE]
  indlist = list()
  for (i in 1:nrow(indmat))
    indlist[[i]] = indmat[i,]
  if ("multicore" %in% search()) { #(is.loaded("mc_fork", PACKAGE="multicore")) {
  mclapply( indlist, function(indvec) {
    i = indvec[1]
    j = indvec[2]
      cat("sms", i, "chr", j, "\n")
      for (k in 1:ngenes) {
       ex <<- exprs(listOfSms[[i]])[k,]
       gfmlaList[[i]][[2]] = as.name("ex")
       tmpc = snp.rhs.tests( gfmlaList[[i]], snp.data=smList(listOfSms[[i]])[[j]], 
          data=pData(listOfSms[[i]]), family="gaussian", ...)@chisq
       if (fillNA) {
          isna = which(is.na(tmpc))
          if (length(isna)>0) 
             tmpc[isna] = rchisq(length(isna), 1)
          }
       if (vmode != "short") shortfac = 1.0
       fflist[[j]][,k,add=TRUE] = tmpc*shortfac
       }  # end k
      }, mc.cores=ncores, mc.set.seed=mc.set.seed)  # end j/mclapply
   }
   else {
    lapply( indlist, function(indvec) {
      i = indvec[1]
      j = indvec[2]
        cat("sms", i, "chr", j, "\n")
        for (k in 1:ngenes) {
         ex <<- exprs(listOfSms[[i]])[k,]
         gfmlaList[[i]][[2]] = as.name("ex")
         tmpc = snp.rhs.tests( gfmlaList[[i]], snp.data=smList(listOfSms[[i]])[[j]], 
            data=pData(listOfSms[[i]]), family="gaussian", ...)@chisq
         if (fillNA) {
            isna = which(is.na(tmpc))
            if (length(isna)>0) 
               tmpc[isna] = rchisq(length(isna), 1)
            }
         if (vmode != "short") shortfac = 1.0
         fflist[[j]][,k,add=TRUE] = tmpc*shortfac
         }  # end k
        })  # end j/mclapply
   } # end windows adaptation
   names(fflist) = chrnames
   ans = list(fflist=fflist, call=theCall, runname=runname, targdir=targdir, 
     generangetag=generangetag,
     filenames=filenames, df=nsms, 
     vmode=vmode, shortfac=shortfac, sessionInfo=sessionInfo(), wd=getwd(), expdataList=expd)
   assign(runname, new("multffManager", ans))
   save(list=runname, file=paste(runname, ".rda", sep=""))
   invisible(get(runname))
}
  

diagffCC = function (sms, gfmla, targdir = ".", runname = "foo", overwriteFF = TRUE, 
    ncores = 2, vmode = "short", shortfac = 100, mc.set.seed=TRUE, fillNA=TRUE, ...) 
{
  if (!file.exists(targdir)) stop("targdir must exist prior to invocation of multffCT")
    theCall = match.call()
    if (!is(sms, "smlSet")) 
        stop("need smlSet as first arg")
    fnhead = paste(targdir, "/", runname, "_", sep = "")
    expd = experimentData(sms)
    nchr = length(smList(sms))
    nsnps = sapply(smList(sms), ncol)
    chrnames = names(smList(sms))
    require(annotation(sms), character.only = TRUE)
    annlib = gsub(".db", "", annotation(sms))
    rmap = revmap(get(paste(annlib, "CHR", sep = "")))
    diaglist = list()
    generangetags = list()
    for (i in chrnames) {
        gs = intersect(featureNames(sms), get(as.character(i), rmap))
        diaglist[[i]] = sms[chrnum(i), ]
        diaglist[[i]] = sms[probeId(gs), ]
        generangetags[[i]] = paste(gs[1], "_", gs[length(gs)], 
            sep = "")
    }
    genenamelist = lapply(diaglist, function(x) featureNames(x))
    ngenelist = lapply(genenamelist, length)
    rslist = lapply(smList(sms), colnames)
    filenames = lapply(chrnames, function(x) paste(fnhead, "snpsOnChr", 
        x, "_", generangetags[[x]], "_df", 1, ".ff", sep = ""))
    fflist = lapply(1:nchr, function(x) ff(initdata = 0, dim = c(nsnps[x], 
        ngenelist[[x]]), dimnames = list(rslist[[x]], genenamelist[[x]]), 
        vmode = vmode, filename = filenames[[x]], overwrite = overwriteFF))
    if ("multicore" %in% search()) { #(is.loaded("mc_fork", PACKAGE="multicore")) {
     mclapply(1:nchr, function(curchr) {
        cat("chr", curchr, "\n")
        cursms = diaglist[[curchr]]
        ngenes = length(featureNames(cursms))
        for (k in 1:ngenelist[[curchr]]) {
            ex <<- exprs(cursms)[k, ]
            gfmla[[2]] = as.name("ex")
            tmpc = snp.rhs.tests(gfmla, snp.data = smList(cursms)[[curchr]],
                data = pData(cursms), family = "gaussian", ...)@chisq
            if (fillNA) {
                isna = which(is.na(tmpc))
                if (length(isna) > 0) 
                  tmpc[isna] = rchisq(length(isna), 1)
            }
            if (vmode != "short") 
                shortfac = 1
            fflist[[curchr]][, k, add = TRUE] = tmpc * shortfac
        }
    }, mc.cores = ncores, mc.set.seed = mc.set.seed)
   }
   else {
    lapply(1:nchr, function(curchr) {
        cat("chr", curchr, "\n")
        cursms = diaglist[[curchr]]
        ngenes = length(featureNames(cursms))
        for (k in 1:ngenelist[[curchr]]) {
            ex <<- exprs(cursms)[k, ]
            gfmla[[2]] = as.name("ex")
            tmpc = snp.rhs.tests(gfmla, snp.data = smList(cursms)[[curchr]],
                data = pData(cursms), family = "gaussian", ...)@chisq
            if (fillNA) {
                isna = which(is.na(tmpc))
                if (length(isna) > 0) 
                  tmpc[isna] = rchisq(length(isna), 1)
            }
            if (vmode != "short") 
                shortfac = 1
            fflist[[curchr]][, k, add = TRUE] = tmpc * shortfac
        }
      })
    }  # end windows adaptation
    names(fflist) = chrnames
    ans = list(fflist = fflist, call = theCall, runname = runname, 
        targdir = targdir, generangetags = generangetags, filenames = filenames, 
        df = 1, vmode = vmode, shortfac = shortfac, sessionInfo = sessionInfo(), 
        wd = getwd(), expdataList = expd)
    assign(runname, new("multffManager", ans))
    save(list = runname, file = paste(runname, ".rda", sep = ""))
    invisible(get(runname))
}
