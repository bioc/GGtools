
setGeneric("oneFit", function(racExSet, geneid, snpid, fitfun, ...)
  standardGeneric("oneFit"))

getpsid = function( geneid, annostring, one.only=TRUE ) {
  require(annostring, character.only=TRUE)
  gnlist = as.list(get(paste(annostring,"SYMBOL",sep="")))
  fonly = sapply(gnlist, function(x)x[1]) # sometimes there are multiple symbols for aprobeset
  find = which(fonly == geneid)
  if (one.only) {
     psid = names(fonly)[find[1]]
     if (length(find) > 1) warning(paste("multiple probesets match", geneid, 
	"using", psid))
     }
  else psid = names(fonly)[find]
  psid
}

setClass("genesym", contains="character")
genesym = function(x) new("genesym", x)

setMethod("oneFit", c("racExSet", "genesym", "character", "function"),
  function(racExSet, geneid, snpid, fitfun, ...) {
    psid = getpsid( geneid, annotation(racExSet) )
    snpvals = snps(racExSet)[snpid,]
    exvals = exprs(racExSet)[psid,]
    ndf = data.frame(exvals,snpvals)
    names(ndf) = c(geneid, snpid)
    fmla = as.formula(paste(as.character(geneid), "~",
		as.character(snpid)))
    fitfun(fmla, data=ndf)
})

setMethod("oneFit", c("racExSet", "character", "character", "function"),
  function(racExSet, geneid, snpid, fitfun, ...) {
    snpvals = snps(racExSet)[snpid,]
    exvals = exprs(racExSet)[geneid,]
    ndf = data.frame(exvals,snpvals)
    names(ndf) = c(geneid, snpid)
    fmla = as.formula(paste(as.character(geneid), "~",
		as.character(snpid)))
    fitfun(fmla, data=ndf)
})

setMethod("oneFit", c("racExSet", "list", "character", "function"),
  function(racExSet, geneid, snpid, fitfun, ...) {
    snpvals = snps(racExSet)[snpid,]
    exvals = geneid
    if (length(geneid) > 1) stop("geneid must be list of length 1")
    geneName = names(geneid)
    if (length(geneName) == 0) stop("geneid must have names attr")
    exvals = geneid[[1]]
    if (!is.numeric(exvals)) stop("geneid list content must be numeric")
    geneid = geneName
## as before below
    if (length(exvals) != length(snpvals)) stop("lengths of numeric expression data and snp allele counts do not agree")
    ndf = data.frame(exvals,snpvals)
    names(ndf) = c(geneid, snpid)
    fmla = as.formula(paste(as.character(geneid), "~",
		as.character(snpid)))
    fitfun(fmla, data=ndf)
})

setMethod("oneFit", c("racExSet", "list", "formula", "function"),
  function(racExSet, geneid, snpid, fitfun, ...) {
#    fixFmla = function(x) {
#      x = as.character(x)
#      x = x[x != "~"]
#      x = strsplit(x, "\\+")
#      extractVars = function(x) gsub(" ","",gsub(")","", gsub("factor\\(","",x)))
#      as.character(sapply(x, extractVars))
#      }
    if (length(as.list(snpid))>2) stop("formula must be of form ~..., no dep. var.")
    snpvars = all.vars(snpid)
    snpvals = snps(racExSet)[snpvars,]
    if (length(snpvars)>1) snpvals=t(snpvals)
    exvals = geneid
    if (length(geneid) > 1) stop("geneid must be list of length 1")
    geneName = names(geneid)
    if (length(geneName) == 0) stop("geneid must have names attr")
    exvals = geneid[[1]]
    if (!is.numeric(exvals)) stop("geneid list content must be numeric")
    geneid = geneName
## as before below
#    if (length(exvals) != length(snpvals)) stop("lengths of numeric expression data and snp allele counts do not agree")
    ndf = data.frame(exvals,snpvals)
    names(ndf) = c(geneid, snpvars)
# now build the formula required by fitfun out of the formula fragment and depvar name
    snpid[[3]] = snpid[[2]]
    snpid[[2]] = as.name(substitute(geneid))
    fitfun(snpid, data=ndf)
})
