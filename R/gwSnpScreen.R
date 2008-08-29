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
 wrapex = function(x) paste("probeId(\"", x, "\")", sep="")
 toks = geneIds(gs)
 idty = geneIdType(gs)
 if (is(idty, "SymbolIdentifier"))
   resps = wrapg(toks)
 else if (is(idty, "AnnotationIdentifier"))
   resps = wrapex(toks)
 else stop("geneIdType for gene set must be either SymbolIdentifier or AnnotationIdentifier")
 lapply(paste(resps, pred, sep="~"), formula)
}


setGeneric("gwSnpTests", function( sym, sms, cnum, ...) standardGeneric("gwSnpTests"))

setMethod("gwSnpTests", c("formula", "smlSet", "cnumOrMissing"),
  function( sym, sms, cnum, ...) {
    if (!missing(cnum)) {
      if (length(cnum) != 1) stop("only supports scalar chrnum cnum at present")
      sms = sms[cnum,]
    }
    theCall = match.call()
#
# decode formula -- note that for a gene set we are essentially recursing
# so it is hard to factor this part out
#
    respObj = eval(sym[[2]]) # we know sym is a formula, sym[[2]] is dep var
    if (is(respObj, "genesym")) {
      annpack = sms@annotation
      require(annpack, character.only=TRUE)
      rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
      pid = AnnotationDbi::get( as(respObj, "character"), rmap )
      if (length(pid) == 0) stop(paste("cannot map", respObj, "in", annpack, sep=""))
      pid = intersect(pid, featureNames(sms))
      if (length(pid) > 1) {
        warning(paste("several probes/sets map to", respObj, "; using", pid[1], sep=""))
        print(pid) 
        pid = pid[1]
        }
      }
    else if (is(respObj, "probeId")) pid = respObj
    else if (is(respObj, "GeneSet")) {
       fms = gsetFmla2FmlaList(sym)
       theCall = match.call()
       if (!missing(cnum)) ans = lapply(fms, function(z) gwSnpTests(z, sms, cnum, ...))
       else ans = lapply(fms, function(z) gwSnpTests(z, sms, ...))
       names(ans) = geneIds(respObj)
       return(new("multiGwSnpScreenResult", geneset=respObj, call=theCall, ans))
       }
    else stop("response in formula must be of class genesym, probeId, or GeneSet")
#
# at this point we have the featureName that we need
#
    pname = as.character(respObj)
    assign(pname, exprs(sms)[pid,]) # expression phenotype genename
    alld = data.frame(get(pname), pData(sms))
    names(alld)[1] = pname
    sym[[2]] = as.name(pname)  # replace the dependent variable spec in fmla
    allsst = lapply( smList(sms), function(x) snp.rhs.tests(sym, family="gaussian",
        snp.data=x, data=alld))
    testType = "Gaussian"
# as of 8 july, we have data frames instead of snp.tests.single objects
# need to coerce
    mksts = function(x) { 
      new("snp.tests.single", chisq=cbind(`1 df`=x$Chi.squared, `2 df`=NA),
          snp.names=rownames(x), N=x$Df.residual+x$Df, N.r2=numeric(0))
    }
    allsst = lapply(allsst, mksts)
#
# return cwSnpScreenResult if chromosome specific, otherwise gwSnpScreenResult
#
    if (!missing(cnum)) return(new("cwSnpScreenResult", gene=respObj, psid=pid,
         annotation=sms@annotation, chrnum=cnum, 
	 call=theCall, testType= testType, allsst))
    new("gwSnpScreenResult", gene=respObj, psid=pid,
         annotation=sms@annotation, 
         call=theCall, testType= testType, allsst)
    })


setMethod("gwSnpTests", c("formula", "smlSet", "snpdepth"),
  function( sym, sms, cnum, ...) {
    if (cnum < 250) stop("with snpdepth numeric third argument you are defining the number of best snps to save per chromosome; it must exceed 250\n")
    theCall = match.call()
    respObj = eval(sym[[2]]) # we know sym is a formula, sym[[2]] is dep var
    if (is(respObj, "genesym")) {
      annpack = sms@annotation
      require(annpack, character.only=TRUE)
      rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
      pid = get( as(respObj, "character"), rmap )
      if (length(pid) == 0) stop(paste("cannot map", respObj, "in", annpack, sep=""))
      pid = intersect(pid, featureNames(sms))
      if (length(pid) > 1) {
        warning(paste("several probes/sets map to", respObj, "; using", pid[1], sep=""))
        print(pid) 
        pid = pid[1]
        }
      }
    else if (is(respObj, "probeId")) pid = respObj
    else if (is(respObj, "GeneSet")) {
       fms = gsetFmla2FmlaList(sym)
       theCall = match.call()
       ans = lapply(fms, function(z) gwSnpTests(z, sms, ...))
       names(ans) = geneIds(respObj)
       tmp <- new("multiGwSnpScreenResult", geneset=respObj, call=theCall, ans)
       names(tmp@.Data) = geneIds(respObj)
       return( filterSnpTests( tmp, cnum ) )
       }
    else if (is(respObj, "chrnum")) {
       rmap = revmap(get(paste(sms@annotation, "CHR", sep="")))
       allpid = get(respObj, rmap)
       allsym = unlist(mget(allpid, get(paste(sms@annotation, "SYMBOL", sep=""))))
	stop("not done")
#       fms = gsetFmla2FmlaList(sym)
#       theCall = match.call()
#       ans = lapply(fms, function(z) gwSnpTests(z, sms, ...))
#       names(ans) = geneIds(respObj)
#       tmp <- new("multiGwSnpScreenResult", geneset=respObj, call=theCall, ans)
#       names(tmp@.Data) = geneIds(respObj)
#       return( filterSnpTests( tmp, cnum ) )
       }
    else stop("response in formula must be of class genesym, probeId, or GeneSet")
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
    tmp = new("gwSnpScreenResult", gene=respObj, psid=pid,
         annotation=sms@annotation, call=theCall,
	 allsst)
    return( filterSnpTests( tmp, cnum ) )
    })

gwSnpScreen = function(...) {
 .Deprecated("gwSnpTests", package="GGtools", "use gwSnpTests instead of gwSnpScreen")
}
