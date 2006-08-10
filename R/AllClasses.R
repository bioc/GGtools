
setClass("ggExprSet", contains="exprSet")

 setMethod("show", "ggExprSet",
   function(object ) {
      dm <-dim(exprs(object))
      ngenes <- dm[1]
      nsamples <- dm[2]
      cat("GG Expression Set (exprSet catering for many SNP attributes) with \n\t", ngenes, " genes\n\t", sep="")
      cat(nsamples, "samples\n")
      cat("There are ", np <- ncol(pData(object)), " attributes; names include:\n")
      cat(colnames(pData(object))[1:min(c(5,np))],"\n")
   }
)

setGeneric("snps", function(x) standardGeneric("snps"))
setMethod("snps", "ggExprSet", function(x) pData(phenoData(x)))

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

setClass("gtexSet", contains="eSet")
make_gtexSet = function(exprs, snps, pd, mi, anno) {
    if (!is(exprs, "matrix")) 
        stop("exprs must be of class matrix")
    if (!is(snps, "matrix")) 
        stop("snps must be of class matrix")
    if (!is(pd, "phenoData")) 
        stop("pd must be of class phenoData")
    aa = assayDataNew("lockedEnvironment", exprs=exprs, snps=snps)
    new("gtexSet", assayData = aa, phenoData = pd, experimentData = mi, 
        annotation = anno)
}

