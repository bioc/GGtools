# we want to inherit all available infrastructure in a combination
# of ExpressionSet and SnpCallSet.  We have an arbitrary number of SnpCallSets to
# associate
setClass("oGtypeExSet", 
  representation(calls="list"), contains="ExpressionSet")

# the major problem is to harmonize the phenoData objects associated with all
# the components.  We assume that the combine() method of Biobase will
# flag important inconsistencies.  We mangle up varLabels so that we know
# where they came from (append an sty or nsp...)
buildOGTES = function(exSet, ...) {
 if (!(is(exSet, "ExpressionSet"))) stop("first argument must be ExpressionSet instance")
 extra = list(...)
 enames = names(extra)
 if (is.null(enames)) stop("extra args to buildOGTES must have names")
 ncallSets = length(extra)
 cls = sapply(extra, function(x) is(x, "SnpCallSet"))
 if (!(all(cls))) stop(
   paste("all extra elements must inherit from SnpCallSet, but got", paste(cls, ,collapse=" ")))
 pds = lapply(extra, phenoData)
 annos = sapply(extra, annotation)
 tags = gsub(".*\\.", "", annos)
 epd = phenoData(exSet)
 nex = length(pds)
 for (i in 1:nex) 
   varLabels(pds[[i]]) = paste(varLabels(pds[[i]]), tags[i], sep=".")
 annos = c(annotation(exSet), annos)
 names(annos) = c("exprs", annos[-1])
 fullPD = phenoData(exSet)
 for (i in 1:nex) 
   fullPD = combine(fullPD, pds[[i]])
 list(fullPD, annos)
 new("oGtypeExSet", exprs=exprs(exSet), phenoData=fullPD, annotation=annos, calls=extra)
}
 
setMethod("show", "oGtypeExSet", function(object) {
 cat("instance of oGtypeExSet (combination of expression\n    data and oligo-based SnpCallSets)\n")
# from eSet show:
  adim <- dim(object)
  if (length(adim)>1)
  cat("assayData for expression:",
      if (length(adim)>1) paste(adim[[1]], "features,", adim[[2]], "samples") else NULL,
      "\n")
  cat("  element names:", paste(assayDataElementNames(object), collapse=", "), "\n")
 cat("assayData for SNPs is collected from", nc <- length(object@calls), "SnpCallSet instances.\n")
 nn = names(object@calls)
 for (i in 1:length(nn)) {
   cat("SnpCallSet", nn[i], ":\n")
   cat(paste("    ", dim(object@calls[[i]])[1], "features.\n"))
   }
 cat("phenoData combined over expression and SNP contributions:\n")
 show(phenoData(object))
})
