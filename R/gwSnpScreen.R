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


setGeneric("gwSnpTests", function( sym, sms, cnum, cs, ...) standardGeneric("gwSnpTests"))

setMethod("gwSnpTests", c("formula", "smlSet", "cnumOrMissing", "missing"),
  function( sym, sms, cnum, cs, ...) {
    if (!missing(cnum)) {
      if (length(cnum) != 1) stop("only supports scalar chrnum cnum at present")
      sms = sms[cnum,]
    }
    theCall = match.call()
    infmla = sym
#
# decode formula -- note that for a gene set we are essentially recursing
# so it is hard to factor this part out
#
    respObj = eval(sym[[2]]) # we know sym is a formula, sym[[2]] is dep var
    if (is(respObj, "phenoVar")) { pid = as.character(respObj@.Data) }
    else if (is(respObj, "genesym")) {
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
       if (!missing(cnum)) ans = lapply(fms, function(z) {
               if (options()$verbose) cat(".")
               gwSnpTests(z, sms, cnum, ...)
               })
       else ans = lapply(fms, function(z) {
          if (options()$verbose) cat(".")
          gwSnpTests(z, sms, ...)
          })
       if (options()$verbose) cat("\n")
       names(ans) = geneIds(respObj)
       SI = new("SessionInfo", sessionInfo())
       if (!missing(cnum)) ans = lapply(ans, function(x) { x@chrnum = cnum; x })
       return(new("multiGwSnpScreenResult", geneset=respObj, call=theCall, 
          sessionInfo = SI, ans))
       }
    else stop("response in formula must be of class phenoVar, genesym, probeId, or GeneSet")
#
# at this point we have the featureName that we need
#
    pname = as.character(respObj)
    sym[[2]] = as.name(pname)  # replace the dependent variable spec in fmla
    if (!is(respObj, "phenoVar")) {
      assign(pname, exprs(sms)[pid,]) # expression phenotype genename
      alld = data.frame(get(pname), pData(sms))
      names(alld)[1] = pname
      }
    else alld = pData(sms)
    allsst = lapply( smList(sms), function(x) snp.rhs.tests(sym, family="gaussian",
        snp.data=x, data=alld))
    testType = "Gaussian"
## as of 8 july, we have data frames instead of snp.tests.single objects
## need to coerce
#    mksts = function(x) { 
#      new("snp.tests.single", chisq=cbind(`1 df`=x$Chi.squared, `2 df`=NA),
#          snp.names=rownames(x), N=x$Df.residual+x$Df, N.r2=numeric(0))
#    }
#    allsst = lapply(allsst, mksts)
#
# return cwSnpScreenResult if chromosome specific, otherwise gwSnpScreenResult
#
    SI = new("SessionInfo", sessionInfo())
    if (!missing(cnum)) return(new("cwSnpScreenResult", gene=respObj, psid=pid,
         annotation=sms@annotation, chrnum=cnum, 
	 call=theCall, sessionInfo=SI, testType= testType, allsst)) # modFmla=infmla, allsst))
    new("gwSnpScreenResult", gene=respObj, psid=pid,
         annotation=sms@annotation, sessionInfo=SI,
         call=theCall, testType= testType, allsst)
    })


setMethod("gwSnpTests", c("formula", "smlSet", "snpdepth", "missing"),
  function( sym, sms, cnum, cs, ...) {
    if (cnum < 250) stop("with snpdepth numeric third argument you are defining the number of best snps to save per chromosome; it must exceed 250\n")
    theCall = match.call()
    respObj = eval(sym[[2]]) # we know sym is a formula, sym[[2]] is dep var
    if (is(respObj, "phenoVar")) { pid = as.character(respObj@.Data)     }
    else if (is(respObj, "genesym")) {
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
       SI = new("SessionInfo", sessionInfo())
       tmp <- new("multiGwSnpScreenResult", geneset=respObj, call=theCall, 
            sessionInfo=SI, ans)
       names(tmp@.Data) = geneIds(respObj)
       return( filterSnpTests( tmp, cnum ) )
       }
    else if (is(respObj, "chrnum")) {
#
# in this segment we transform chrnum spec to gene set
# and then reinvoke
#
       require( sms@annotation, character.only=TRUE, quietly=TRUE )
    if (is(respObj, "phenoVar")) { pid = as.character(respObj@.Data) }
       require( GSEABase, quietly=TRUE )
       rmap = revmap(get(paste(gsub(".db", "", sms@annotation), "CHR", sep="")))
       allpid = get(as(respObj,"character"), rmap)
       allsym = unlist(mget(allpid, get(paste(gsub(".db", "", sms@annotation), "SYMBOL", sep=""))))
       gs = GeneSet(unique(allsym), geneIdType=SymbolIdentifier())
       tmp = as.list(sym)
       tmp[[2]] = gs
       sym = as.formula(tmp)
       return( gwSnpTests( sym, sms, cnum, ...) )
       }
    else stop("response in formula must be of class phenoVar, genesym, probeId, or GeneSet")
    pname = as.character(respObj)
    sym[[2]] = as.name(pname)  # replace the dependent variable spec in fmla
    if (!is(respObj, "phenoVar")) {
      assign(pname, exprs(sms)[pid,]) # expression phenotype genename
      alld = data.frame(get(pname), pData(sms))
      names(alld)[1] = pname
      }
    else alld = pData(sms)
    allsst = lapply( smList(sms), function(x) snp.rhs.tests(sym, family="gaussian",
        snp.data=x, data=alld))
# as of 8 july, we have data frames instead of snp.tests.single objects
# need to coerce
#    mksts = function(x) { 
#      new("snp.tests.single", chisq=cbind(`1 df`=x$Chi.squared, `2 df`=NA),
#          snp.names=rownames(x), N=x$Df.residual+x$Df, N.r2=numeric(0))
#    }
#    allsst = lapply(allsst, mksts)
    SI = new("SessionInfo", sessionInfo())
    tmp = new("gwSnpScreenResult", gene=respObj, psid=pid,
         annotation=sms@annotation, call=theCall, sessionInfo=SI,
	 allsst)
    return( filterSnpTests( tmp, cnum ) )
    })

setGeneric("residTests", function(fit, sms, litfmla, rsnum) 
  standardGeneric("residTests"))

setMethod("residTests", c("cwSnpScreenResult", "smlSet", "formula", "missing"), function(fit, sms, litfmla, rsnum) {
  theCall = match.call()
  top = rownames(topSnps(fit))[1]
  smm = smList(sms)[[fit@chrnum]][, top, drop=FALSE]
  baseRAC = as( smm, "numeric" )
  ex = exprs(sms)[ fit@psid, ]
  ok = 1:length(ex)
  bad = NULL
  if (any(lkna <- is.na(baseRAC))) bad = which(lkna)
  if (any(lkna <- is.na(ex))) bad = union(bad, which(lkna))
  if (length(bad)>0) ok = ok[-bad]
  res = resid(lm(ex ~ baseRAC, subset=ok))
  #fmla = fit@formula
  litfmla[[2]] = as.name("res")
  alld = data.frame(res, pData(sms)[ok,])
  allsst = lapply( smList(sms), function(x) snp.rhs.tests(litfmla, family="gaussian",
        snp.data=x, data=alld))
# mksts = function(x) {
##
## I THINK YOU NEED TO DROP
##
#      new("snp.tests.single", chisq=cbind(`1 df`=x$Chi.squared, `2 df`=NA),
#          snp.names=rownames(x), N=x$Df.residual+x$Df, N.r2=numeric(0))
#    }
#    allsst = lapply(allsst, mksts)
#
# return cwSnpScreenResult if chromosome specific, otherwise gwSnpScreenResult
#
    SI = new("SessionInfo", sessionInfo())
    return(new("cwSnpScreenResult", gene=fit@gene, psid=fit@psid,
         annotation=sms@annotation, chrnum=fit@chrnum, sessionInfo=SI,
         call=theCall, testType= "Gaussian resid", allsst)) # modFmla=fit@formula, allsst))

})
 
