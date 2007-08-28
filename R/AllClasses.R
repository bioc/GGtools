# GGtools AllClasses.R (c) 2006 VJ Carey

# helper classes to figure out semantics of character strings

setClass("snpID", contains="character")
snpID = function(x) new("snpID", x)

setClass("genesym", contains="character")
genesym = function(x) new("genesym", x)

setGeneric("snps", function(x) standardGeneric("snps"))

# helper class for snp screen output

setClass("snpScreenResult", representation(call="call", gene="character", 
    locs="numeric", chr="character", fitter="ANY"), contains="list")

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
