# GGtools infrastructure from hm2ceu, Mar 31 2008 (c) VJ Carey

valsml = function(object) {
 allns = sapply(smList(object), nrow)
 if (!(all(allns==allns[1]))) 
    return("varying numbers of rows in elements of smList")
 if ((sl <- length(smList(object))) != (cl <- length(object@chromInds)))
    return(paste("length of chromInds vector [", cl, "] not identical to that of smList(object) [",
      sl, "]"))
 nna = names(annotation(object))
 if ((length(nna) != 2) || (!(all(nna == c("exprs", "snps")))))
    return("annotation slot must be vector with names 'exprs' and 'snps'")
# if (!(is(annotation(object)[2], "snpLocNCref")))
#    return("second item in annotation slot must inherit from snpLocNCref")
 return(TRUE)
}

snpLocPathClo = function(pkg="GGdata", subdir="extdata") function(x) {
  dir(system.file(subdir, package=pkg), full=TRUE, patt=x)
}

setClass("smlSet", contains="eSet", 
   representation(smlEnv="environment", snpLocPathMaker="function",
     chromInds="numeric", organism="character", snpLocPackage="character",
     snpLocRef="character", activeSnpInds="numeric"),
   validity=valsml, prototype=prototype(
       new("VersionedBiobase",
               versions=c(classVersion("eSet"), smlSet="1.0.0")),
           phenoData = new("AnnotatedDataFrame",
             data=data.frame(),
             varMetadata=data.frame(
               labelDescription=character(0))),
	   annotation=character(0),
	   smlEnv = {e = new.env(); assign("smList", list(), e); e},
 	   snpLocPathMaker=snpLocPathClo(), snpLocRef=character(0),
   	   snpLocPackage=character(0), activeSnpInds=numeric(0),
           chromInds = numeric(0)))
	   
setClass("gwSnpScreenResult", contains="list",
   representation(gene="character", psid="character", annotation="character",
      snpLocPackage="character", snpLocNCDFref="character",
      activeSnpInds="numeric"))

setClass("multiGwSnpScreenResult", representation(geneset="GeneSet"),
   contains="list")


setClass("chrnum", contains="numeric")
setClass("rsNum", contains="character")

setClassUnion("cnumOrMissing", c("chrnum", "missing"))

setGeneric("chrnum", function(x) standardGeneric("chrnum"))
setMethod("chrnum", "numeric", function(x) new("chrnum", x))
setGeneric("rsNum", function(x) standardGeneric("rsNum"))
setMethod("rsNum", "character", function(x) new("rsNum", x))
setClass("genesym", contains="character")
setGeneric("genesym", function(x) standardGeneric("genesym"))
setMethod("genesym", "character",  function(x) new("genesym", x))

setClass("cwSnpScreenResult", contains="gwSnpScreenResult",
   representation(chrnum="chrnum"))

setClass("snpLocNCref", contains="character")
setMethod("show", "snpLocNCref", function(object) {
 cat(gsub("..*/", "...", object), " [netCDF]\n")
})


# GGtools AllClasses.R (c) 2006 VJ Carey -- legacy below

setClass("snpMeta", 
  representation(meta="environment", chromosome="character"))

setClass("snpMetaWhole",  contains="snpMeta",
  representation(chrbounds="numeric", chrlabs="character"))

# helper classes to figure out semantics of character strings

setClass("snpID", contains="character")
snpID = function(x) new("snpID", x)

setClass("exFeatID", contains="character")
exFeatID = function(x) new("exFeatID", x)

setClass("genesym", contains="character")
genesym = function(x) new("genesym", x)

setGeneric("snps", function(x, chr) standardGeneric("snps"))

# helper class for snp screen output

setClass("snpScreenResult", representation(call="call", gene="character", 
    locs="numeric", chr="character", fitter="ANY",
      annotation="character"), contains="list")

setClass("twSnpScreenResult", representation(call="call", genes="character", 
    locs="numeric", fitter="ANY"), contains="list")

setMethod("show", "twSnpScreenResult", function(object) {
   cat("twSnpScreenResult\n")
   cat("call: ")
   print(object@call)
   cat("Genes (selection):\n", selectSome(names(object)))
   cat("\nFirst fit object:\n")
   cat("---\n")
   show(object[[1]])
   cat("--- [there are", length(object)-1, "more]\n")
   })

topSnpsOLD = function(x, n=10) {
  lapply(x, function(x) sort(extract_p(x))[1:n])
}

# key class for genetical genomics Aug 2006 -- rare allele
# count combined with expression (racExSet)

setClass("racExSet", representation(
    racAssays="AssayData",
    rarebase="character", SNPalleles="character"), #, snpMeta="snpMeta"), 
    contains="eSet",
    prototype = prototype(racAssays=assayDataNew()))

setMethod("initialize", "racExSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   exprs = new("matrix"),
		   racs = new("matrix"),
		   rarebase = character(),
                   SNPalleles = character()) {
            .Object = callNextMethod(.Object,
                           assayData = assayDataNew(
                             exprs=exprs),
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
	    .Object@racAssays = assayDataNew(racs =racs)
            .Object@SNPalleles = SNPalleles
	    .Object@rarebase = rarebase
            .Object
          })

setClass("GGfitter", representation(name="character", func="function"))

fastAGMfitter = new("GGfitter", name="fastAGM", func=fastAGM)
fastHETfitter = new("GGfitter", name="fastHET", func=fastHET)
