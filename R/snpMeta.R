setClass("snpMeta", 
  representation(meta="environment", chromosome="character"))

setClass("snpMetaWhole",  contains="snpMeta",
  representation(chrbounds="numeric", chrlabs="character"))

setGeneric("chromosome", function(x) standardGeneric("chromosome"))
setMethod("chromosome", "snpMeta", function(x) x@chromosome)

df2snpMeta = function(df, chrom) {
 mm = new.env()
 assign("meta", df, mm)
 new("snpMeta", meta=mm, chromosome=chrom)
}

setMethod("[", "snpMeta", function(x, i, j, ..., drop=FALSE) {
 if (missing(drop)) drop=FALSE
 if (missing(j)) get("meta", x@meta)[i,,drop=drop]
 else if (missing(i)) get("meta", x@meta)[,j,drop=drop]
 else get("meta", x@meta)[i,j,drop=drop]
})

setMethod("show", "snpMeta", function(object) {
 cat("snp metadata for chromosome", chromosome(object), "\n")
 cat("first five records:\n")
 print(object[1:5,])
})

setAs("snpMeta", "data.frame", function(from) get("meta", from@meta))
setMethod("dim", "snpMeta", function(x) dim(as(x, "data.frame")))
setMethod("nrow", "snpMeta", function(x) nrow(as(x, "data.frame")))


