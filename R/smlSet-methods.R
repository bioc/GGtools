	   
#setGeneric("smlEnv", function(x) standardGeneric("smlEnv"))
#setMethod("smlEnv", "smlSet", function(x) x@smlEnv)
#setGeneric("smList", function(x) standardGeneric("smList"))
#setMethod("smList", "smlSet", function(x) x@smlEnv$smList)
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
 cat("snp locs package:", object@snpLocPackage, "; SQLite ref:", 
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
      tmp = exprs(x)[,j,drop=FALSE]
      x@assayData=assayDataNew("lockedEnvironment", exprs=tmp)
      }
 if (!missing(i) && is(i, "chrnum")) {
      L2 = L2[match(i, x@chromInds)]
      x@chromInds = x@chromInds[match(i, x@chromInds)]
      }
 else if (!missing(i) && is(i, "exFeatID")) {
      tmp = exprs(x)[i,,drop=FALSE]
      x@assayData=assayDataNew("lockedEnvironment", exprs=tmp)
      x@featureData = x@featureData[i,]
      }
 else if (!missing(i) && is(i, "genesym")) {
      annpack = x@annotation["exprs"]
      require(annpack, character.only=TRUE)
      rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
      psid = get(i, rmap)
      if (length(psid) > 1) warning("gene symbol matches multiple probe sets, using first")
      psid = psid[1]
      if (missing(j)) return(x[ exFeatID(psid), ])
      else return(x[ exFeatID(psid), j])
      }
 else if (!missing(i)) stop("only exFeatID (probe), chrnum (chrom for SNPs),
    or genesym (HUGO)  row-selections for smlSet supported")
 eL = new.env(hash=TRUE)
 assign("smList", L2, eL)
 x@smlEnv = eL
 #if (!missing(j)) phenoData(x)@data = phenoData(x)@data[j,,drop=FALSE]
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

setMethod("gwSnpScreen", c("probeId", "smlSet", "cnumOrMissing"),
  function( sym, sms, cnum, ...) {
    if (!missing(cnum)) {
      if (length(cnum) != 1) stop("only supports scalar chrnum cnum at present")
      sms = sms[cnum,]
    }
    theCall = match.call()
    if (length(sym) > 1) stop("for multiple gene expression analysis, please use a GSEABase::GeneSet instance")
    ph = exprs(sms)[sym,]
    allsst = lapply( smList(sms), function(x) single.snp.tests(pheno=ph,
        snp.data=x))
    pid = sym
    if (!missing(cnum)) return(new("cwSnpScreenResult", gene=sym, psid=pid,
         annotation=sms@annotation, chrnum=cnum, 
         snpLocPackage=sms@snpLocPackage,
	 snpLocExtRef=sms@snpLocRef, activeSnpInds=sms@activeSnpInds, call=theCall,
         allsst))
    new("gwSnpScreenResult", gene=sym, psid=pid,
         annotation=sms@annotation, 
	 snpLocPackage=sms@snpLocPackage, snpLocExtRef=
           sms@snpLocRef, activeSnpInds=sms@activeSnpInds, call=theCall,
	 allsst)
    })


setMethod("gwSnpScreen", c("genesym", "smlSet", "cnumOrMissing"),
  function( sym, sms, cnum, ...) {
    if (!missing(cnum)) {
      if (length(cnum) != 1) stop("only supports scalar chrnum cnum at present")
      sms = sms[cnum,]
    }
    theCall = match.call()
    if (length(sym) > 1) stop("for multiple gene expression analysis, please use a GSEABase::GeneSet instance")
    annpack = sms@annotation["exprs"]
    require(annpack, character.only=TRUE)
    rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
    pid = get( as(sym, "character"), rmap )
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
	 snpLocExtRef=sms@snpLocRef, activeSnpInds=sms@activeSnpInds, call=theCall,
         allsst))
    new("gwSnpScreenResult", gene=sym, psid=pid,
         annotation=sms@annotation, 
	 snpLocPackage=sms@snpLocPackage, snpLocExtRef=
           sms@snpLocRef, activeSnpInds=sms@activeSnpInds, call=theCall,
	 allsst)
    })

setGeneric("getSnpLocs", function(x,filterActive) standardGeneric("getSnpLocs"))
setGeneric("getSnpChroms", function(x,filterActive) standardGeneric("getSnpChroms"))

setMethod("getSnpLocs", c("smlSet", "missing"), function(x, filterActive) {
 getSnpLocs(x, TRUE)
})

setMethod("getSnpLocs", c("SQLiteConnection", "character"), function(x, filterActive) {
 RSQLite:::sqliteQuickColumn(x, filterActive, "loc")
})

setGeneric("getSnpLocConn", function(x)standardGeneric("getSnpLocConn"))
setMethod("getSnpLocConn", "smlSet", function(x) {
  get(x@snpLocRef)()
})

setMethod("getSnpLocs", c("smlSet", "logical"), function(x, filterActive) {
 locs = getSnpLocs( getSnpLocConn(x), gsub("_dbconn", "", x@snpLocRef))
 if (!filterActive) return(locs)
 else {
      activeChr = x@chromInds
      chrnum = getSnpChroms( getSnpLocConn(x), gsub("_dbconn", "", x@snpLocRef) )
      if (isTRUE(all.equal(sort(activeChr), sort(unique(as.numeric(chrnum))))))
            return(locs)
      else return(locs[which(chrnum %in% activeChr)])
      }
})
       
setMethod("getSnpChroms", c("SQLiteConnection", "character"), function(x, filterActive) {
 RSQLite:::sqliteQuickColumn(x, filterActive, "chrnum")
})

setMethod("getSnpChroms", c("smlSet", "missing"), function(x, filterActive) {
 getSnpChroms(x, TRUE)
})

setMethod("getSnpChroms", c("smlSet", "logical"), function(x, filterActive) {
 chrnum = getSnpChroms( getSnpLocConn(x), gsub("_dbconn", "", x@snpLocRef))
 if (!filterActive) return(chrnum)
 else {
      activeChr = x@chromInds
      if (isTRUE(all.equal(sort(activeChr), sort(unique(as.numeric(chrnum))))))
            return(chrnum)
      else return(chrnum[which(chrnum %in% activeChr)])
      }
})


setGeneric("getAlleles", function(x, rs, ...) standardGeneric("getAlleles"))
setMethod("getAlleles", c("smlSet", "rsNum"), function (x, rs)
{
    allrs = lapply(smList(x), colnames)
    #hits = lapply(allrs, function(x) grep(rs, x))
    hits = lapply(allrs, function(z) which(z == rs))
    kpi = sapply(hits, function(z) length(z) > 0)
    if (!(any(kpi))) stop("rs number not found in columns of smlSet")
    ans = list(chr = names(kpi[kpi]), col = hits[[which(kpi)]])
    as(smList(x)[[ans$chr]][, ans$col], "character")
})


setGeneric("plot_EvG", function(gsym, rsn, sms, ...)standardGeneric("plot_EvG"))
setMethod("plot_EvG", c("genesym", "rsNum", "smlSet"), 
    function(gsym, rsn, sms, ...) {
      annpack = sms@annotation["exprs"]
      require(annpack, character.only=TRUE)
      rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
      psid = get(gsym, rmap)
      if (length(psid) > 1) warning("gene symbol matches multiple probe sets, using first")
      psid = psid[1]
      ex = exprs(sms)[psid, ] # this returns a data.frame!?! for hmyriB36
      if (is(ex, "data.frame")) ex = as.numeric(ex)
      gt = factor(getAlleles( sms, rsn ))
      plot(ex~gt, ylab=gsym, xlab=rsn, ...)
      points(jitter(as.numeric(gt),.4), ex, col="gray", pch=19)
      })
setMethod("plot_EvG", c("probeId", "rsNum", "smlSet"), 
    function(gsym, rsn, sms, ...) {
#      annpack = sms@annotation["exprs"]
#      library(annpack, character.only=TRUE)
#      rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
      psid = gsym
      if (length(psid) > 1) warning("gene symbol matches multiple probe sets, using first")
      psid = psid[1]
      ex = exprs(sms)[psid, ] # this returns a data.frame!?! for hmyriB36
      if (is(ex, "data.frame")) ex = as.numeric(ex)
      gt = factor(getAlleles( sms, rsn ))
      plot(ex~gt, ylab=gsym, xlab=rsn, ...)
      points(jitter(as.numeric(gt),.4), ex, col="gray", pch=19)
      })
setMethod("plot_EvG", c("character", "character", "smlSet"), 
    function(gsym, rsn, sms, ...) {
      warning("arguments not cast using, e.g., genesym, assuming HUGO symbol and dbSnp IDs are passed")
      annpack = sms@annotation["exprs"]
      require(annpack, character.only=TRUE)
      rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
      psid = get(gsym, rmap)
      if (length(psid) > 1) warning("gene symbol matches multiple probe sets, using first")
      psid = psid[1]
      ex = exprs(sms)[psid, ] # this returns a data.frame!?! for hmyriB36
      if (is(ex, "data.frame")) ex = as.numeric(ex)
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

#setClass("multiGwSnpScreenResult", representation(geneset="GeneSet"),
#   contains="list")
setMethod("gwSnpScreen", c("GeneSet", "smlSet", "cnumOrMissing"),
  function( sym, sms, cnum, ...) {
  gid = as(geneIds(sym), "character")
  ng = length(gid)
  out = list()
  for (i in 1:ng) {
   if (!missing(cnum)) {
    out[[i]] = try(gwSnpScreen(genesym(gid[i]), sms, cnum, ...))
    if (inherits(out[[i]], "try-error")) out[[i]] =  
               try(gwSnpScreen(probeId(gid[i]), sms, cnum, ...))
    }
   else {
    out[[i]] = try(gwSnpScreen(genesym(gid[i]), sms, ...))
    if (inherits(out[[i]], "try-error")) out[[i]] =  
               try(gwSnpScreen(probeId(gid[i]), sms, ...))
    }
   }
  if (options()$verbose) cat(sym, "\n")
  new("multiGwSnpScreenResult", geneset=sym, out)
})

setMethod("show", "multiGwSnpScreenResult", function(object) {
 cat("multi genome-wide snp screen result:\n")
 cat("gene set used as response:\n")
 show(object@geneset)
 cat("there are", length(object), "results.\n")
 cat("the call was:\n")
 print(object@call)
})

gsetFmla2FmlaList = function(fm) {
#
# this function takes a formula with a GeneSet on the lhs
# and returns a list of gene-specific formulas
#
 flist = as.list(fm)
 gs = eval(flist[[2]])
 if (length(flist[[3]]) > 1) pred = paste(flist[[3]][-1], collapse="+")
 else if (length(flist[[3]]) == 1) pred = flist[[3]]
 if (!is(gs, "GeneSet")) stop("needs GeneSet instance in response position")
 wrapg = function(x) paste("genesym(\"", x, "\")", sep="")
 wrapex = function(x) paste("exFeatID(\"", x, "\")", sep="")
 toks = geneIds(gs)
 idty = geneIdType(gs)
 if (is(idty, "SymbolIdentifier"))
   resps = wrapg(toks)
 else if (is(idty, "AnnotationIdentifier"))
   resps = wrapex(toks)
 else stop("geneIdType for gene set must be either SymbolIdentifier or AnnotationIdentifier")
 lapply(paste(resps, pred, sep="~"), formula)
}


setMethod("gwSnpScreen", c("formula", "smlSet", "cnumOrMissing"),
  function( sym, sms, cnum, ...) {
    if (!missing(cnum)) {
      if (length(cnum) != 1) stop("only supports scalar chrnum cnum at present")
      sms = sms[cnum,]
    }
    theCall = match.call()
    respObj = eval(sym[[2]]) # we know sym is a formula, sym[[2]] is dep var
    if (is(respObj, "genesym")) {
      annpack = sms@annotation["exprs"]
      require(annpack, character.only=TRUE)
      rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
      pid = get( as(respObj, "character"), rmap )
      if (length(pid) == 0) stop(paste("cannot map", respObj, "in", annpack, sep=""))
      if (length(pid) > 1) {
        warning(paste("several probes/sets map to", respObj, "; using", pid[1], sep=""))
        print(pid) 
        pid = pid[1]
        }
      }
    else if (is(respObj, "exFeatID")) pid = respObj
    else if (is(respObj, "GeneSet")) {
       fms = gsetFmla2FmlaList(sym)
       theCall = match.call()
       if (!missing(cnum)) ans = lapply(fms, function(z) gwSnpScreen(z, sms, cnum, ...))
       else ans = lapply(fms, function(z) gwSnpScreen(z, sms, ...))
       return(new("multiGwSnpScreenResult", geneset=respObj, call=theCall, ans))
       }
    else stop("response in formula must be of class genesym, exFeatID, or GeneSet")
    pname = as.character(respObj)
    assign(pname, exprs(sms)[pid,]) # expression phenotype genename
    alld = data.frame(get(pname), pData(sms))
    names(alld)[1] = pname
    sym[[2]] = as.name(pname)  # replace the dependent variable spec in fmla
    allsst = lapply( smList(sms), function(x) snp.rhs.tests(sym, family="gaussian",
        snp.data=x, data=alld))
# as of 8 july, we have data frames instead of snp.tests.single objects
# need to coerce
    mksts = function(x) { 
      new("snp.tests.single", chisq=cbind(`1 df`=x$Chi.squared, `2 df`=NA),
          snp.names=rownames(x), N=x$Df.residual+x$Df, N.r2=numeric(0))
    }
    allsst = lapply(allsst, mksts)
    if (!missing(cnum)) return(new("cwSnpScreenResult", gene=as.character(respObj), psid=pid,
         annotation=sms@annotation, chrnum=cnum, 
         snpLocPackage=sms@snpLocPackage,
	 snpLocExtRef=sms@snpLocRef, activeSnpInds=sms@activeSnpInds, call=theCall,
         allsst))
    new("gwSnpScreenResult", gene=as.character(respObj), psid=pid,
         annotation=sms@annotation, 
	 snpLocPackage=sms@snpLocPackage, snpLocExtRef=
           sms@snpLocRef, activeSnpInds=sms@activeSnpInds, call=theCall,
	 allsst)
    })
