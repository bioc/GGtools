
#setClass("ggExprSet", contains="exprSet")
#
# setMethod("show", "ggExprSet",
#   function(object ) {
#      dm <-dim(exprs(object))
#      ngenes <- dm[1]
#      nsamples <- dm[2]
#      cat("GG Expression Set (exprSet catering for many SNP attributes) with \n\t", ngenes, " genes\n\t", sep="")
#      cat(nsamples, "samples\n")
#      cat("There are ", np <- ncol(pData(object)), " attributes; names include:\n")
#      cat(colnames(pData(object))[1:min(c(5,np))],"\n")
#   }
#)

setGeneric("snps", function(x) standardGeneric("snps"))
#setMethod("snps", "ggExprSet", function(x) pData(phenoData(x)))

#> getClass("eSet")
#Virtual Class
#
#Slots:
#                                                               
#Name:           assayData          phenoData     experimentData
#Class:          AssayData AnnotatedDataFrame              MIAME
#                                            
#Name:          annotation  .__classVersion__
#Class:          character           Versions
#
#Extends: 
#Class "VersionedBiobase", directly
#Class "Versioned", by class "VersionedBiobase"
#
#Known Subclasses: "ExpressionSet", "MultiSet", "SnpSet"
#> getClass("AssayData")
#Extended class definition ( "ClassUnionRepresentation" )
#Virtual Class
#
#No Slots, prototype of class "NULL"
#
#Known Subclasses: "list", "environment"
#


#setValidity("ggAssayData", function(object) {
# msg <- NULL
# nn <- names(object)
# if (!all.equal(c("exprs", "snps"), nn)) msg <-
#     "ggAssayData needs exprs and snps only"
# if (is.null(msg)) return(TRUE)
# else return(msg)
#})

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


setMethod("snps", "racExSet", function(x) get("racs",x@racAssays))
setMethod("exprs", "racExSet", function(object) get("exprs",object@assayData))
setGeneric("racAssays", function(x) standardGeneric("racAssays"))
setMethod("racAssays", "racExSet", function(x) x@racAssays)
setGeneric("snpNames", function(x) standardGeneric("snpNames"))
setMethod("snpNames", "racExSet", function(x) featureNames(x@racAssays))
setGeneric("rarebase", function(x) standardGeneric("rarebase"))
setMethod("rarebase", "racExSet", function(x) x@rarebase)
setGeneric("SNPalleles", function(x) standardGeneric("SNPalleles"))
setMethod("SNPalleles", "racExSet", function(x) x@SNPalleles)

setMethod("show", "racExSet", function(object) {
    cat("racExSet instance (SNP rare allele count + expression)\n")
    cat("rare allele count assayData:\n")
  cat("  Storage mode:", storageMode(racAssays(object)), "\n")
  nms <- selectSome(snpNames(object))
  cat("  featureNames:", paste(nms, collapse=", "))
  if ((len <- length(snpNames(object))) > length(nms))
    cat(" (", len, " total)", sep="")
  cat("\n  Dimensions:\n")
  print(Biobase:::assayDataDims(racAssays(object)))
  cat("\nexpression assayData\n")
  cat("  Storage mode:", storageMode(object), "\n")
  nms <- selectSome(featureNames(object))
  cat("  featureNames:", paste(nms, collapse=", "))
  if ((len <- length(featureNames(object))) > length(nms))
    cat(" (", len, " total)", sep="")
  cat("\n  Dimensions:\n")
  print(dims(object))
  cat("\nphenoData\n")
  show(phenoData(object)) 
  cat("\n")
  show(experimentData(object))
  cat("\nAnnotation ")
  show(annotation(object))
    })

make_racExSet = function(exprs, racs, rarebase, SNPalleles, pd, mi, anno) {
    if (!is(exprs, "matrix")) 
        stop("exprs must be of class matrix")
    if (!is(racs, "matrix")) 
        stop("racs must be of class matrix")
    if (!is(pd, "phenoData")) 
        stop("pd must be of class phenoData")
    names(SNPalleles) = rownames(racs)
    names(rarebase) = rownames(racs)
    new("racExSet", exprs=exprs, racs=racs, rarebase=rarebase, 
        SNPalleles = SNPalleles,
        phenoData = pd, experimentData = mi, 
        annotation = anno)
}

