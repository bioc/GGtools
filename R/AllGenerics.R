setGeneric("racAssays", function(x) standardGeneric("racAssays"))
setGeneric("snpNames", function(x) standardGeneric("snpNames"))
setGeneric("rarebase", function(x) standardGeneric("rarebase"))
setGeneric("SNPalleles", function(x) standardGeneric("SNPalleles"))

setGeneric("snpScreen", function(racExSet, snpMeta, gene, formTemplate, fitter, gran, ...)
    standardGeneric("snpScreen"))
