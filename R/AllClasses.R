# GGtools AllClasses.R (c) 2006 VJ Carey

# helper classes to figure out semantics of character strings

setClass("snpID", contains="character")
snpID = function(x) new("snpID", x)

setClass("exFeatID", contains="character")
exFeatID = function(x) new("exFeatID", x)

setClass("genesym", contains="character")
genesym = function(x) new("genesym", x)

setGeneric("snps", function(x) standardGeneric("snps"))

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

topSnps = function(x, n=10) {
  lapply(x, function(x) sort(extract_p(x))[1:n])
}

# key class for genetical genomics Aug 2006 -- rare allele
# count combined with expression (racExSet)

setClass("racExSet", representation(
    racAssays="AssayData",
    rarebase="character", SNPalleles="character"), contains="eSet",
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
