#setGeneric("filterSnpTests", function(x, n) standardGeneric(
# "filterSnpTests"))
#
#setMethod("filterSnpTests",
#  "multiGwSnpScreenResult", function(x, n) {
#    tmp = lapply( x@.Data, filterGWS, n )
#    gc()
#    names(tmp) = GSEABase::geneIds(x@geneset)
#    new("filteredMultiGwSnpScreenResult", geneset=x@geneset,
#      call=x@call, tmp)
#})
#
#setMethod("filterSnpTests",
#  "gwSnpScreenResult", function(x, n) {
#   filterGWS(x, n)
#})
#
#
#setMethod("plot", "filteredGwSnpScreenResult", function(x, y, ...) {
# pp = lapply(x@.Data, p.value) 
# boxplot(lapply(pp, function(x)-log10(x)), main=x@gene, xlab="chromosome",
#   ylab="-log10 p [GLM]")
# xx = try(require(org.Hs.eg.db, quietly=TRUE))
# if (!inherits(xx, "try-error")) {
#  if (is(x@gene, "genesym")) rmap = revmap(org.Hs.egSYMBOL)
# else if (is(x@gene, "probeId"))  {
#           require(x@annotation, character.only=TRUE, quietly=TRUE)
#           rmap = get(paste(gsub(".db", "", x@annotation), "ENTREZID", sep=""))
#           }
#
#  else {
#    warning("x@gene is neither symbol nor probeID, we do not plot the location.")
#    return(invisible(NULL))
#    }
#    egid = get(x@gene, rmap)
#    ch = try(get(egid, org.Hs.egCHR))
#    if (!inherits(ch, "try-error")) {
#       if (ch == "X") ch = 23
#       else if (ch == "Y") ch = 24
#       axis(3, at=as.numeric(ch), col="red", labels=" ")
#       }
#    }
#})

#setMethod("plot", "filteredMultiGwSnpScreenResult", function(x, y, ...) {
# stop("please select the desired gene-specific result via [[ and plot directly\n")
#})
geneRanges = function(ids, annopkg, extend=0) {
 require(annopkg, character.only=TRUE)
 anbase = gsub(".db", "", annopkg)
 chrmap = get(paste(anbase, "CHR", sep=""))
 stmap = get(paste(anbase, "CHRLOC", sep=""))
 endmap = get(paste(anbase, "CHRLOCEND", sep=""))
 chrs = mget(ids, chrmap)
 chrs = unlist(sapply(chrs, "[", 1))
 chrs = paste("chr", chrs, sep="")
 sts = mget(ids, stmap)
 sts = unlist(sapply(sts, "[", 1))
 ens = mget(ids, endmap)
 ens = unlist(sapply(ens, "[", 1))
 negsts = which(sts < 0)
 sts[negsts] = -sts[negsts]
 ens[negsts] = -ens[negsts]
 st = pmax(1,sts-extend)
 en = ens+extend
 st[is.na(sts)] = 1
 en[is.na(sts)] = 2
 RangedData(IRanges(st,en), space=chrs, name=ids)
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
       names(ans) = GSEABase::geneIds(respObj)
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
        snp.data=x, data=alld, uncertain=TRUE))
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
       names(ans) = GSEABase::geneIds(respObj)
       SI = new("SessionInfo", sessionInfo())
       tmp <- new("multiGwSnpScreenResult", geneset=respObj, call=theCall, 
            sessionInfo=SI, ans)
       names(tmp@.Data) = GSEABase::geneIds(respObj)
       return( filterSnpTests( tmp, cnum ) )
       }
    else if (is(respObj, "chrnum")) {
#
# in this segment we transform chrnum spec to gene set
# and then reinvoke
#
       require( sms@annotation, character.only=TRUE, quietly=TRUE )
    if (is(respObj, "phenoVar")) { pid = as.character(respObj@.Data) }
       #require( GSEABase, quietly=TRUE )
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
    allsst = lapply( smList(sms), function(x) snp.rhs.tests(sym, family="gaussian", uncertain=TRUE,
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
        snp.data=x, data=alld, uncertain=TRUE))
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
 



setClass("chunksize", contains="numeric")
chunksize = function(x) new("chunksize", as.numeric(x))


setMethod("gwSnpTests", c("formula", "smlSet", "snpdepth", "chunksize"),
 function(sym, sms, cnum, cs) {
# assumes a gene set is response of formula
  theCall = match.call(call=sys.call(2))
  gn = GSEABase::geneIds(gs <- eval(sym[[2]]))
  ng = length(gn)
   chunklabs = function (n, chunksize) 
   {
       bas = 1:n
       tool = ceiling(n/chunksize)
       as.numeric(cut(bas, tool))
   }
  gspl = split(gn, chunklabs(ng, cs))
  csets = lapply( gspl, function(x) gs[x] )
  savesym = sym
  out = list()
  for (i in 1:length(gspl)) {
      nsym = savesym
      nsym[[2]] = csets[[i]]
      out[[i]] = gwSnpTests(nsym, sms, cnum)
      gc()
      }
# this list has all the tests filtered already, so filterSnpTests is not 
# needed
  flattened = unlist(out, recursive=FALSE)
  names(flattened) = gn
  ans = new("filteredMultiGwSnpScreenResult", geneset=gs,
      call=theCall, flattened)
  names(ans@.Data) = gn
  ans
  })

#setClass("maxchisq", contains="list")
#setMethod("show", "maxchisq", function(object) {
# cat("GGtools maxchisq structure.\n")
# cat("The call was:\n")
# print(object$theCall)
# cat("The original call in multffManager was:\n")
# print(object$mgrcall)
# cat("Excerpt:\n")
# print(lapply(object[c("maxchisq", "bestFeats")], function(x) head(x[[1]])))
#})

   
#setGeneric("min_p_vals", function(mcs, mtcorr, type, sidedness) standardGeneric("min_p_vals"))
#setMethod("min_p_vals", c("maxchisq", "character", "character", "numeric"), function(mcs, mtcorr, type, sidedness) {
# sidedness = as.integer(sidedness)
# if (sidedness != 1 & sidedness != 2) stop("sidedness must be 1 or 2")
# pv = lapply(mcs$maxchisq, function(x) pmin(1, sidedness*(1-pchisq(x, mcs$df))))
# npv = lapply(mcs$maxchisq, names)
# mtcorrp = function(x, mtcorr) {
#   tmp = mt.rawp2adjp(x, mtcorr)
#   tmp$adjp[ order(tmp$index), mtcorr ]
# }
# if (mtcorr == "none") adjpv = pv
# else {
#   stop("owing to a namespace complication, please use mtcorr = 'none' and compute corrections on your own.")
#   #require(multtest)
#   if (type == "chr_specific")
#     adjpv = lapply( pv, function(x) mtcorrp(x, mtcorr))
#   else if (type=="global") {
#     ulp = unlist(pv)
#     uln = unlist(npv)
#     names(ulp) = uln
#     adjpv = mtcorrp(ulp, mtcorr)
#     names(adjpv) = uln
#     anslist = list()
#     for (i in 1:length(npv)) {
#       anslist[[i]] = adjpv[ npv[[i]] ] # restore chromosomal list structure
#       names(anslist[[i]]) = npv[[i]]
#     }
#     names(anslist) = names(npv)
#     return(anslist)
#   }  
# }
# for (i in 1:length(adjpv)) names(adjpv[[i]]) = npv[[i]]
# adjpv
#})

#setClass("multffManager", contains="list")
#setMethod("show", "multffManager", function(object) {
# require(ff, quietly=TRUE)
# cat("multffManager instance. The call was:\n")
# print(object$call)
# cat("There are ", length(object$filenames), " ff files.\n")
# cat("Excerpt from first file:\n")
# ngenes = ncol(object$fflist[[1]])
# print(object$fflist[[1]][1:4,1:min(4,ngenes)])
#})

#setMethod("[", c("multffManager", "rsid", "missing"), function(x, i, j, ..., drop=TRUE) {
#  div = 1.0
#  if (x$vmode == "short") div = x$shortfac
#  rsn = rsnum(x)
#  pres = which(sapply(rsn, function(y) any(i %in% y)))
#  if (length(pres) < 1) stop("rsid not found in rownames of any ff matrix from multffCT")
#  ffl = x$fflist[pres]
#  tmp = lapply(ffl, function(y)y[intersect(rownames(y),as(i,"character")),,drop=FALSE])
#  lapply(tmp, function(z)z/div)
#})

#setMethod("[", c("multffManager", "missing", "probeId"), function(x, i, j, ..., drop=TRUE) {
#  div = x$shortfac
#  if (x$vmode != "short") div = 1.0
#  if (length(j) > 50 & !isTRUE(getOption("ggt_manyGenes"))) stop("to request more than 50 genes with [ please set option ggt_manyGenes to TRUE with options(ggt_manyGenes=TRUE)\nand recognize that you may be creating a large RAM image")
#  lapply(x$fflist, function(x)x[, as(j, "character"), drop=FALSE]/div)
#})
#
#setMethod("[", c("multffManager", "rsid", "probeId"), function(x, i, j, ..., drop=TRUE) {
#  tmp = x[ i, , drop=FALSE ]
#  lapply(tmp, function(z) z[, j, drop=FALSE])
#})


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


chkeman = function(object){
# eqtlTestsManager validity test
 allgn = lapply(fflist(object), colnames)
 n1 = allgn[[1]]
 chk = sapply(allgn[-1], function(x)all.equal(x,n1))
 if (!all(chk)) return("fflist colnames not common to all elements")
 if (is.null(names(fflist(object)))) return("fflist elements lack names")
 return(TRUE)
}

# elements of a multffManager list
#> names(dem)
# [1] "fflist"       "call"         "runname"      "targdir"      "generangetag"
# [6] "filenames"    "df"           "vmode"        "shortfac"     "sessionInfo" 
#[11] "wd"           "expdataList" 


setClass("eqtlTestsManager",
 representation(fflist="list", call="call", sess="ANY",
	exdate="ANY", shortfac="numeric", geneanno="character", df="numeric",
        summaryList="list"),
        validity=chkeman)

setAs("multffManager", "eqtlTestsManager", function(from) {
 new("eqtlTestsManager", fflist=from$fflist, call=from$call,
      sess=from$sessionInfo, shortfac=from$shortfac, df=from$df,
      exdate=paste("converted:", date()), geneanno="please supply")
})

setGeneric("probeNames", function(x) standardGeneric("probeNames"))
setMethod("probeNames", "eqtlTestsManager", function(x) {
 colnames(fflist(x)[[1]])  # assumes common gene set for fflist
})


setGeneric("shortfac", function(x)standardGeneric("shortfac"))
setMethod("shortfac", "eqtlTestsManager", function(x)
  x@shortfac)
setGeneric("fflist", function(x)standardGeneric("fflist"))
setMethod("fflist", "eqtlTestsManager", function(x)
  x@fflist)
setGeneric("exdate", function(x)standardGeneric("exdate"))
setMethod("exdate", "eqtlTestsManager", function(x)
  x@exdate)

setMethod("show", "eqtlTestsManager", function(object) {
 cat("eqtlTools results manager, computed", exdate(object), "\n")
 cat("gene annotation:", object@geneanno, "\n")
 cat("There are", length(fflist(object)), "chromosomes analyzed.\n")
 cat("some genes (out of ", length(colnames(fflist(object)[[1]])),"): ", paste(selectSome(colnames(fflist(object)[[1]])),collapse=" "), "\n", sep="")
 cat("some snps (out of ", sum(sapply(fflist(object),nrow)),  "): ", paste(selectSome(rownames(fflist(object)[[1]])),collapse=" "), "\n", sep="")
})


setMethod("[", c("eqtlTestsManager", "rsid", "probeId"),
 function(x, i, j, ..., drop=FALSE) {
#
# ultimately this may not be exposed, serving only for deep
# testing, because a director database may be required for every
# manager
#
 m1 = snpIdMap( as(i, "character"), x )
 ans = lapply(1:length(m1), function(i) fflist(x)[[names(m1)[i]]][ m1[[i]], 
    as(j, "character"), drop=FALSE]/shortfac(x))
 names(ans) = names(m1)
 ans
})

setMethod("[", c("eqtlTestsManager", "missing", "probeId"),
 function(x, i, j, ..., drop=FALSE) {
#
#
 ll = length(fflist(x))
 ans = lapply(1:ll, function(i) fflist(x)[[i]][ , 
    as(j, "character"), drop=FALSE]/shortfac(x))
 names(ans) = names(fflist(x))
 ans
})

setMethod("[", c("eqtlTestsManager", "rsid", "missing"),
 function(x, i, j, ..., drop=FALSE) {
 m1 = snpIdMap( as(i, "character"), x )
 ans = lapply(1:length(m1), function(i) fflist(x)[[names(m1)[i]]][ m1[[i]], 
    , drop=FALSE]/shortfac(x))
 names(ans) = names(m1)
 ans
})

setMethod("[", c("eqtlTestsManager", "ANY", "ANY"),
 function(x, i, j, ..., drop=FALSE) {
 stop("[ for eqtlTestsManager only defined for signature ('rsid', 'probeId') [one may be omitted]")
 })

# director for group of managers

chkmgrs = function(object) {
   mcl = sapply(mgrs(object), class)
   chkc = sapply(mgrs(object), function(x) is(x, "eqtlTestsManager"))
   if (!all(chkc)) return("mgrs slot must only contain list of entities inheriting from eqtlTestsManager")
   annos = sapply(mgrs(object), function(x)x@geneanno)
   if (!all(annos==annos[1])) return("managers do not have identical gene annotation source")
   sids = lapply(mgrs(object), snpIdList)
   slchk = sapply(sids, function(x) all.equal(x, sids[[1]]))
   if (!all(sapply(slchk,isTRUE))) return("managers do not have identical SNP lists")
   return(TRUE)
}

 
   
setClass("cisTransDirector", 
  representation(mgrs="list", indexdbname="character", 
   shortfac="numeric", snptabname="character", probetabname="character", probeanno="character", snptabref="ANY", probetabref="ANY"),
   validity=chkmgrs)

setGeneric("mgrs", function(x) standardGeneric("mgrs"))
setMethod("mgrs", "cisTransDirector", function(x) x@mgrs)
setGeneric("nsnps", function(x) standardGeneric("nsnps"))
setMethod("nsnps", "cisTransDirector", function(x) sum(sapply(fflist(mgrs(x)[[1]]), nrow)))
setGeneric("ngenes", function(x) standardGeneric("ngenes"))
setMethod("ngenes", "cisTransDirector", function(x) sum(sapply(mgrs(x), function(y) ncol(fflist(y)[[1]]))))

nsnp = function(cd) sum(sapply(cd@mgrs[[1]]@fflist, nrow))# function(x) sum(sapply(x@fflist, nrow)))
#ngenes = function(cd) sum(sapply(cd@mgrs, function(x)ncol(x@fflist[[1]])) )


setMethod("show", "cisTransDirector", function(object) {
 cat("eqtlTools cisTransDirector instance.\n")
 cat("there are", length(mgrs(object)), "managers.\n")
 cat("Total number of SNP: ", nsnp(object), "; total number of genes: ", ngenes(object), "\n")
 cat("First:\n")
 show(mgrs(object)[[1]])
 cat("---\n")
 cat("use [ (rsnumvec), (geneidvec) ] to obtain chisq stats; topFeats(), etc.\n")
})


#setMethod("[", c("cisTransDirector", "character", "character"),
#  function(x, i, j, ..., drop=FALSE) {
##    if (length(j)>1) stop("currently only handle single probe reference")
#    snpListChr = unique(as.character(x@snptabref[i,]))
#    if (length(snpListChr)>1) stop("currently only collecting scores for SNP on a single chromosome")
#    probeListEl = sort(unique(as.integer(x@probetabref[j,])))
##
## following assumes common SNP over managers
##
#    mgrlist = lapply(probeListEl, function(z) mgrs(x)[[ z ]])
##    names(mgrlist) = j
#    ans = lapply(1:length(mgrlist), 
#       function(z) fflist(mgrlist[[z]])[[snpListChr]][ i, j[z] ]/shortfac(mgrlist[[z]]))
##    names(ans) = j
#    ans
#})
#


setMethod("[", c("cisTransDirector", "character", "character"),
 function (x, i, j, ..., drop = FALSE) 
 {
    # following will be an index, so numeric
    snpListChr = unique(x@snptabref[i, ])
    if (length(snpListChr) > 1) 
        stop("currently only collecting scores for SNP on a single chromosome")
    prinds = as.integer(x@probetabref[j, ])
    spids = split(j, prinds)
    probeListEl = sort(unique(as.integer(x@probetabref[j, ])))
    if (!all.equal(as.integer(names(spids)), probeListEl)) 
		stop("split of gene names by director element indices has unexpected result")
    mgrlist = lapply(probeListEl, function(z) GGtools:::mgrs(x)[[z]])
    applier = lapply
    if ("multicore" %in% search()) applier = mclapply
    ans = applier(1:length(mgrlist), function(z) GGtools:::fflist(mgrlist[[z]])[[snpListChr]][i, 
        spids[[z]],drop=FALSE]/GGtools:::shortfac(mgrlist[[z]]))
    if (length(ans) == 1) return(ans[[1]])
    bans = ans[[1]]
    for (i in 2:length(ans)) bans = cbind(bans, ans[[i]])
    bans
 }
)


#setMethod("[", c("cisTransDirector", "character", "missing"),
#  function(x, i, j, ..., drop=FALSE) {
#    snpListChr = unique(as.character(x@snptabref[i,]))
#    if (length(snpListChr)>1) stop("currently only collecting scores for SNP on a single chromosome")
##
## following assumes common SNP over managers
##
#    mgrlist = mgrs(x)
#    ans = lapply(1:length(mgrlist), 
#       function(z) fflist(mgrlist[[z]])[[snpListChr]][ i, ]/shortfac(mgrlist[[z]]))
#    allsn = rownames(ans[[1]])
#    allgn = unlist(lapply(ans, colnames))
#    nans = t(sapply(ans, function(x)x))
#    rownames(nans) = allsn
#    colnames(nans) = allgn
#    nans
#})

setMethod("[", c("cisTransDirector", "character", "missing"),
  function(x, i, j, ..., drop=FALSE) {
  j = unlist(probeNames(x))
  callGeneric(x, i, j, drop=drop)
})

setMethod("probeNames", "cisTransDirector", function(x) {
 lapply(mgrs(x), probeNames) 
})

setMethod("[", c("cisTransDirector", "missing", "character"),
  function(x, i, j, ..., drop=FALSE) {
    probeListEl = sort(unique(as.integer(x@probetabref[j,])))
#
# following assumes common SNP over managers
#
    mgrlist = lapply(probeListEl, function(z) mgrs(x)[[ z ]])
    nsnps = length(fflist(mgrlist[[1]]))  # nchr??
#    names(mgrlist) = j
    applier = lapply
    if ("multicore" %in% search()) applier = mclapply
    ans = applier(1:length(mgrlist), 
       function(z) lapply(1:nsnps, function(w) fflist(mgrlist[[z]])[[w]][ , j[z] ]/shortfac(mgrlist[[z]])))
#    names(ans) = j
    ans
})

setGeneric("topSnps", function(x, ...) standardGeneric("topSnps"))
setMethod("topSnps", "cwSnpScreenResult", function(x, n=10) {
   pp = p.value(x@.Data[[1]])
   sn = x@.Data[[1]]@snp.names  # no accessor...
   opp = order(pp, decreasing=FALSE) 
   spp = pp[ opp ]
   df = data.frame(p.val=spp)
   rownames(df) = sn[ opp ]
   df[1:n,,drop=FALSE]
})

setMethod("topSnps", "gwSnpScreenResult", function(x, n=10) {
  ts.df = function (w, n = 10) {
   pp = p.value(w)
   sn = w@snp.names  # no accessor...  # don't need list access here
   opp = order(pp, decreasing=FALSE)
   spp = pp[ opp ]
   df = data.frame(p.val=spp)
   rownames(df) = sn[ opp ]
   df[1:n,,drop=FALSE]
   }
  lapply(x, ts.df, n=n)
})

setAs("cwSnpScreenResult", "RangedData", function(from) {
  allp = p.value(from@.Data[[1]]) 
  rs = from@.Data[[1]]@snp.names
  locstr = snpLocs.Hs(chrnum(from@chrnum), rsid(rs))
  loc = locstr["loc",]
  locrs = paste("rs", locstr["rsid",], sep="")
  allp = allp[locrs]
  
  require(org.Hs.eg.db, quietly=TRUE)
  rmap = revmap(org.Hs.egSYMBOL)
  ch = paste("chr", from@chrnum, sep="")

  rd = RangedData(IRanges(loc, loc), type = "snpeff", group = "gws",
    score = as.numeric(-log10(allp)), space = ch, universe = "hg18")
  allp = rd$score # order probably has changed
  bad = is.na(allp) | !is.finite(allp)
  if (any(bad))
    rd = rd[!bad,]
  rd
})

setAs("cwSnpScreenResult", "GRanges", function(from) {
  allp = p.value(from@.Data[[1]]) # , 1) # assume 1df -- must improve
  rs = from@.Data[[1]]@snp.names
  locstr = snpLocs.Hs(chrnum(from@chrnum), rsid(rs))
  loc = locstr["loc",]
  locrs = paste("rs", locstr["rsid",], sep="")
  allp = allp[locrs]
  
  require(org.Hs.eg.db, quietly=TRUE)
  rmap = revmap(org.Hs.egSYMBOL)
  ch = paste("chr", from@chrnum, sep="")

  rd = GRanges(IRanges(loc, loc), type = "snpeff", group = "gws",
    score = as.numeric(-log10(allp)), seqnames = ch, universe = "hg18")
  names(rd) = names(allp)
  allp = elementMetadata(rd)$score # order probably has changed
  bad = is.na(allp) | !is.finite(allp)
  if (any(bad))
    rd = rd[!bad,]
  rd
})

#setMethod("annotation", "eqtlTestsManager", function(x, ...) {
# x@geneanno
#})

setClass("transManager", representation(base="list"))

setMethod("show", "transManager", function(object){
 require(ff, quietly=TRUE)
 basel = object@base
 cat("transManager instance, created", basel$date, "\n", sep=" ")
 cat("dimension of scores component:\n")
 cat(" number of loci checked: ", nrow(basel$scores), 
   "; genes retained: ", ncol(basel$scores), "\n", sep="")
 cat("the call was:\n")
 print(basel$call)
})

.probesManaged = function(mgr,ffind=1) {
 colnames(mgr@fflist[[ffind]])
}

.snpsManaged = function(mgr,ffind=1) {
 rownames(mgr@fflist[[ffind]])
}

setGeneric("probesManaged", function(mgr, ffind)
 standardGeneric("probesManaged"))
setGeneric("snpsManaged", function(mgr, ffind)
 standardGeneric("snpsManaged"))

setMethod("probesManaged", c("eqtlTestsManager",
     "numeric"), function(mgr, ffind=1)
       {
       .probesManaged(mgr, ffind)
       })

setMethod("snpsManaged", c("eqtlTestsManager",
     "numeric"), function(mgr, ffind=1)
       {
       .snpsManaged(mgr, ffind)
       })
