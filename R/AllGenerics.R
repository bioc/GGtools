setGeneric("racAssays", function(x) standardGeneric("racAssays"))
#setGeneric("snpNames", function(x) standardGeneric("snpNames"))
setGeneric("rarebase", function(x) standardGeneric("rarebase"))
setGeneric("SNPalleles", function(x) standardGeneric("SNPalleles"))

setGeneric("snpScreen", function(racExSet, snpMeta, gene, formTemplate, fitter, gran)
    standardGeneric("snpScreen"))

setGeneric("twSnpScreen", function(racExSet, snpMeta, formTemplate, fitter)
    standardGeneric("twSnpScreen"))

setGeneric("snps", function(x, chr) standardGeneric("snps"))

