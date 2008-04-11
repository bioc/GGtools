	   
setGeneric("smlEnv", function(x) standardGeneric("smlEnv"))
setMethod("smlEnv", "smlSet", function(x) x@smlEnv)
setGeneric("smList", function(x) standardGeneric("smList"))
setMethod("smList", "smlSet", function(x) x@smlEnv$smList)
setGeneric("nsamp", function(x) standardGeneric("nsamp"))
setMethod("nsamp", "smlSet", function(x) nrow(smList(x)[[1]]))
setMethod("sampleNames", "smlSet", function(object) rownames(smList(object)[[1]]))

setMethod("show", "smlSet", function(object) {
 cat("snp.matrix-based genotype set:\n")
 cat("number of samples: ", nsamp(object), "\n")
 cat("number of snp.matrix: ", length(smList(object)), "\n")
 cat("annotation:\n")
 cat(" exprs: ")
 cat( object@annotation["exprs"])
 cat("\n snps: ")
# cat(object@annotation["snps"], "[arg to @snpLocPathMaker]\n")
 cat("snp locs package:", object@snpLocPackage, "; ncdf ref:", 
     object@snpLocRef, "\n")
 if (length(dd <- dim(object@assayData$exprs))>0) {
  cat("Expression data:", dd[1], "x", dd[2], "\n") 
 }
 cat("Phenodata: "); show(phenoData(object))
})


setMethod("exprs", "smlSet", function(object) object@assayData$exprs)
setMethod("[", "smlSet", function(x, i, j, ..., drop=FALSE) {
# i indexes reporters, j samples
 if (missing(i) & missing(j)) return(x)
 L2 = smList(x)
 if (!missing(j)) {
      L2 = lapply(L2, function(z) z[j, ,drop=FALSE])
      phenoData(x)@data = phenoData(x)@data[j,,drop=FALSE]
      }
 if (!missing(i) && is(i, "chrnum")) {
      L2 = L2[match(i, x@chromInds)]
      x@chromInds = x@chromInds[match(i, x@chromInds)]
      }
 else if (!missing(i)) stop("only chrnum 'row' selections for snp.matrix lists supported")
 eL = new.env(hash=TRUE)
 assign("smList", L2, eL)
 x@smlEnv = eL
 if (!missing(j)) phenoData(x)@data = phenoData(x)@data[j,,drop=FALSE]
 #callNextMethod() # for now, adequate to deal with phenoData
 x
})


setGeneric("rawSNP", function(x, chrind)standardGeneric("rawSNP"))
setMethod("rawSNP", c("smlSet", "numeric"), 
    function(x, chrind) {
     if (!(chrind %in% x@chromInds)) stop("chrind must be an element of chromInds slot")
     smList(x)[[match(chrind, x@chromInds)]]@.Data
     })

setGeneric("gwSnpScreen", function(sym, sms, cnum, ...)
      standardGeneric("gwSnpScreen"))

#setMethod("gwSnpScreen", c("genesym", "smlSet", "cnumOrMissing"),
#  function( sym, sms, cnum, ...) {
#    annpack = annotation(sms)["exprs"]
#    library(annpack, character.only=TRUE)
#    rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
#    pid = get( sym, rmap )
#    if (length(pid) == 0) stop(paste("cannot map", sym, "in", annpack, sep=""))
#    if (length(pid) > 1) {
#        warning(paste("several probes/sets map to", sym, "; using", pid[1], sep=""))
#        print(pid) 
#        pid = pid[1]
#        }
#    ph = exprs(sms)[pid,]
#    allsst = lapply( smList(sms), function(x) single.snp.tests(pheno=ph,
#        snp.data=x))
#    new("gwSnpScreenResult", gene=sym, psid=pid,
#         annotation=sms@annotation, allsst)
#    })

setMethod("gwSnpScreen", c("genesym", "smlSet", "cnumOrMissing"),
  function( sym, sms, cnum, ...) {
    if (!missing(cnum)) {
      if (length(cnum) != 1) stop("only supports scalar chrnum cnum at present")
      sms = sms[cnum,]
    }
    annpack = sms@annotation["exprs"]
    library(annpack, character.only=TRUE)
    rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
    pid = get( sym, rmap )
    if (length(pid) == 0) stop(paste("cannot map", sym, "in", annpack, sep=""))
    if (length(pid) > 1) {
        warning(paste("several probes/sets map to", sym, "; using", pid[1], sep=""))
        print(pid) 
        pid = pid[1]
        }
    ph = exprs(sms)[pid,]
    allsst = lapply( smList(sms), function(x) single.snp.tests(pheno=ph,
        snp.data=x))
    if (!missing(cnum)) return(new("cwSnpScreenResult", gene=sym, psid=pid,
         annotation=sms@annotation, chrnum=cnum, 
         snpLocPackage=sms@snpLocPackage,
	 snpLocNCDFref=sms@snpLocRef,
         allsst))
    new("gwSnpScreenResult", gene=sym, psid=pid,
         annotation=sms@annotation, 
	 snpLocPackage=sms@snpLocPackage, snpLocNCDFref=
           sms@snpLocRef,
	 allsst)
    })

setGeneric("getSnpLocs", function(x,filterActive) standardGeneric("getSnpLocs"))

setMethod("getSnpLocs", c("smlSet", "missing"), function(x, filterActive) {
 getSnpLocs(x, TRUE)
})

setMethod("getSnpLocs", c("smlSet", "logical"), function(x, filterActive=FALSE) {
#
# it is now assumed that we have exported the ncdf reference
# in object named x@snpLocRef, in package named x@snpLocPackage
#
  require(x@snpLocPackage, character.only=TRUE)
  ref = get(x@snpLocRef, paste("package:", x@snpLocPackage, sep=""))
# if a smlSet has been subset by chromosome for SNP
# then its chromInds are the 'active' chromosomes
# we return the locations corresponding to these
# [ if it has not been subset, then all chromosomes are 'active' ]
 activeChr = x@chromInds
 availchr = get.var.ncdf(ref, "chr")
 loc = get.var.ncdf(ref, "cumloc")
 if (!filterActive) return(loc)
 if (isTRUE(all.equal(sort(activeChr), sort(unique(as.numeric(availchr))))))
    return(loc)
 else return(loc[which(availchr %in% activeChr)])
 })

#setMethod("getSnpLocs", "smlSet", function(x) {
##
## if a smlSet has been subset by chromosome for SNP
## then its chromInds are the 'active' chromosomes
## we return the locations corresponding to these
## [ if it has not been subset, then all chromosomes are 'active' ]
##
# ncpath = x@snpLocPathMaker(x@annotation["snps"])
# oo = open.ncdf( ncpath )
# on.exit(close(oo))
# activeChr = x@chromInds
# annochr = get.var.ncdf(oo, "chr")
# kpinds = which(annochr %in% activeChr)
# annoloc = get.var.ncdf(oo, "cumloc")
# annoloc[kpinds]
#})

#setGeneric("getSnpChroms", function(x) standardGeneric("getSnpChroms"))
#setMethod("getSnpChroms", "smlSet", function(x) {
##
## if a smlSet has been subset by chromosome for SNP
## then its chromInds are the 'active' chromosomes
## we return the chromosome nums corresponding to these
## [ if it has not been subset, then all chromosomes are 'active' ]
##
# ncpath = x@snpLocPathMaker(x@annotation["snps"])
# oo = open.ncdf( ncpath )
# on.exit(close(oo))
# activeChr = x@chromInds
# annochr = get.var.ncdf(oo, "chr")
# kpinds = which(annochr %in% activeChr)
# annochr[kpinds]
#})


setGeneric("getAlleles", function(x, rs, ...) standardGeneric("getAlleles"))
setMethod("getAlleles", c("smlSet", "rsNum"), function (x, rs)
{
    allrs = lapply(smList(x), colnames)
    hits = lapply(allrs, function(x) grep(rs, x))
    kpi = sapply(hits, function(x) length(x) > 0)
    ans = list(chr = names(kpi[kpi]), col = hits[[which(kpi)]])
    as(smList(x)[[ans$chr]][, ans$col], "character")
})


setGeneric("plot_EvG", function(gsym, rsn, sms, ...)standardGeneric("plot_EvG"))
setMethod("plot_EvG", c("genesym", "rsNum", "smlSet"), 
    function(gsym, rsn, sms, ...) {
      annpack = sms@annotation["exprs"]
      library(annpack, character.only=TRUE)
      rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
      psid = get(gsym, rmap)
      if (length(psid) > 1) warning("gene symbol matches multiple probe sets, using first")
      psid = psid[1]
      ex = exprs(sms)[psid, ]
      gt = factor(getAlleles( sms, rsn ))
      plot(ex~gt, ylab=gsym, xlab=rsn, ...)
      points(jitter(as.numeric(gt),.4), ex, col="gray", pch=19)
      })

setMethod("snps", c("smlSet", "chrnum"), function(x, chr) {
  if (length(chr) != 1) stop("chr must have length 1")
  tmp = as(smList(x[chr,])[[1]], "character")
  rownames(tmp) = sampleNames(x)
  colnames(tmp) = colnames(smList(x[chr,])[[1]])
  t(tmp)
})

