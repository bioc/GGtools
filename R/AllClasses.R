# GGtools AllClasses.R (c) 2006 VJ Carey

# helper classes to figure out semantics of character strings

setClass("snpID", contains="character")
snpID = function(x) new("snpID", x)

setClass("genesym", contains="character")
genesym = function(x) new("genesym", x)

setGeneric("snps", function(x) standardGeneric("snps"))

# helper class for snp screen output

setClass("snpScreenResult", representation(call="call", gene="character", 
    locs="numeric", chr="character", fittertok="character"), contains="list")

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

#
#setMethod("snps", "racExSet", function(x) get("racs",x@racAssays))
#setMethod("exprs", "racExSet", function(object) get("exprs",object@assayData))
#
#setGeneric("racAssays", function(x) standardGeneric("racAssays"))
#setMethod("racAssays", "racExSet", function(x) x@racAssays)
#setGeneric("snpNames", function(x) standardGeneric("snpNames"))
#setMethod("snpNames", "racExSet", function(x) featureNames(x@racAssays))
#setGeneric("rarebase", function(x) standardGeneric("rarebase"))
#setMethod("rarebase", "racExSet", function(x) x@rarebase)
#setGeneric("SNPalleles", function(x) standardGeneric("SNPalleles"))
#setMethod("SNPalleles", "racExSet", function(x) x@SNPalleles)
#
#setMethod("show", "racExSet", function(object) {
#    cat("racExSet instance (SNP rare allele count + expression)\n")
#    cat("rare allele count assayData:\n")
#  cat("  Storage mode:", storageMode(racAssays(object)), "\n")
#  nms <- selectSome(snpNames(object))
#  cat("  featureNames:", paste(nms, collapse=", "))
#  if ((len <- length(snpNames(object))) > length(nms))
#    cat(" (", len, " total)", sep="")
#  cat("\n  Dimensions:\n")
#  print(Biobase:::assayDataDims(racAssays(object)))
#  cat("\nexpression assayData\n")
#  cat("  Storage mode:", storageMode(object), "\n")
#  nms <- selectSome(featureNames(object))
#  cat("  featureNames:", paste(nms, collapse=", "))
#  if ((len <- length(featureNames(object))) > length(nms))
#    cat(" (", len, " total)", sep="")
#  cat("\n  Dimensions:\n")
#  print(dims(object))
#  cat("\nphenoData\n")
#  show(phenoData(object)) 
#  cat("\n")
#  show(experimentData(object))
#  cat("\nAnnotation ")
#  show(annotation(object))
#    })
#
#make_racExSet = function(exprs, racs, rarebase, SNPalleles, pd, mi, anno) {
#    if (!is(exprs, "matrix")) 
#        stop("exprs must be of class matrix")
#    if (!is(racs, "matrix")) 
#        stop("racs must be of class matrix")
#    if (!is(pd, "phenoData")) 
#        stop("pd must be of class phenoData")
#    names(SNPalleles) = rownames(racs)
#    names(rarebase) = rownames(racs)
#    new("racExSet", exprs=exprs, racs=racs, rarebase=rarebase, 
#        SNPalleles = SNPalleles,
#        phenoData = pd, experimentData = mi, 
#        annotation = anno)
#}
#
#
