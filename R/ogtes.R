# we want to inherit all available infrastructure in a combination
# of ExpressionSet and SnpCallSet.  We have an arbitrary number of SnpCallSets to
# associate
setClass("oGtypeExSet", 
  representation(snpCalls="list", dbConns="list"), contains="ExpressionSet")

setAs("oGtypeExSet", "ExpressionSet", function(from) {
 anno = as.character(annotation(from)["exprs"])
 new("ExpressionSet", assayData=assayData(from),
      phenoData=phenoData(from), featureData=featureData(from),
      experimentData=experimentData(from), annotation=anno)
})
 


setMethod("initialize", "oGtypeExSet", 
  function(.Object, snpCalls=new("list"), dbConns=new("list"), x) { # x is ExpressionSet
#
# first transfer the ExpressionSet data
#
	.Object = callNextMethod(.Object,
		assayData=assayData(x),
		phenoData=phenoData(x),
		featureData=featureData(x),
		experimentData=experimentData(x),
		annotation=annotation(x)) # what if slots of ExpressionSet change?
#
# now combine all the phenoDatas from exprs and snpCalls
#
	 pds = lapply(snpCalls, phenoData)
	 annos = sapply(snpCalls, annotation)
	 tags = gsub(".*\\.", "", annos)
	 epd = phenoData(x)
	 nex = length(pds)
	 for (i in 1:nex)
	   varLabels(pds[[i]]) = paste(varLabels(pds[[i]]), tags[i], sep=".")
	 annos = c(annotation(x), annos)
	 names(annos) = c("exprs", annos[-1])
	 fullPD = phenoData(x)
	 for (i in 1:nex)
	   fullPD = combine(fullPD, pds[[i]])
#
# finally set up the extended slot and put in the enhanced
# phenoData and annotation
#
	.Object@snpCalls = snpCalls
	.Object@dbConns = dbConns
	.Object@phenoData = fullPD
	.Object@annotation = annos
        .Object
})

# the major problem is to harmonize the phenoData objects associated with all
# the components.  We assume that the combine() method of Biobase will
# flag important inconsistencies.  We mangle up varLabels so that we know
# where they came from (append an sty or nsp...)
#buildOGTES = function(exSet, ...) {
# if (!(is(exSet, "ExpressionSet"))) stop("first argument must be ExpressionSet instance")
# extra = list(...)
# enames = names(extra)
# if (is.null(enames)) stop("extra args to buildOGTES must have names")
# ncallSets = length(extra)
# cls = sapply(extra, function(x) is(x, "SnpCallSet"))
# if (!(all(cls))) stop(
#   paste("all extra elements must inherit from SnpCallSet, but got", paste(cls, ,collapse=" ")))
# pds = lapply(extra, phenoData)
# annos = sapply(extra, annotation)
# tags = gsub(".*\\.", "", annos)
# epd = phenoData(exSet)
# nex = length(pds)
# for (i in 1:nex) 
#   varLabels(pds[[i]]) = paste(varLabels(pds[[i]]), tags[i], sep=".")
# annos = c(annotation(exSet), annos)
# names(annos) = c("exprs", annos[-1])
# fullPD = phenoData(exSet)
# for (i in 1:nex) 
#   fullPD = combine(fullPD, pds[[i]])
## list(fullPD, annos)
# tmp = new("ExpressionSet", assayData=assayData(exSet),
#		phenoData=fullPD,
#		featureData=featureData(exSet),
#		experimentData=experimentData(exSet),
#		annotation=annos)
# new("oGtypeExSet", snpCalls=extra, tmp)
#}
 
setMethod("show", "oGtypeExSet", function(object) {
 cat("instance of oGtypeExSet (combination of expression\n    data and oligo-based SnpCallSets)\n")
# from eSet show:
  adim <- dim(object)
  if (length(adim)>1)
  cat("assayData for expression:",
      if (length(adim)>1) paste(adim[[1]], "features,", adim[[2]], "samples") else NULL,
      "\n")
  cat("  element names:", paste(assayDataElementNames(object), collapse=", "), "\n")
 cat("assayData for SNPs is collected from", nc <- length(object@snpCalls), "SnpCallSet instances.\n")
 nn = names(object@snpCalls)
 if (length(nn) > 0) {
   for (i in 1:length(nn)) {
     cat("SnpCallSet", nn[i], ":\n")
     cat(paste("    ", dim(object@snpCalls[[i]])[1], "features.\n"))
     }
   }
   cat("phenoData combined over expression and SNP contributions:\n")
   show(phenoData(object))
})
