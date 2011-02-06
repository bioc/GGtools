
setClass("egtSet", representation(gtlist="list", snpNames="list"), contains="eSet",
 prototype=prototype(new("ExpressionSet"), gtlist=list()))

convertSS = function(smlSet, ffstub="ffsnp_") {
 ss = smList(smlSet)
 sn = lapply(ss, colnames)
 cn = names(ss)
 fn = paste(ffstub, cn, ".ff", sep="")
 ss = lapply(ss, function(x) {dimnames(x) <- NULL; x})
 ffrefs = lapply(1:length(ss), function(i) ff(ss[[i]], dim=dim(ss[[i]]), vmode="raw",
		overwrite=TRUE, filename=fn[i]))
 tmp = lapply(ffrefs, function(x) close.ff(x))
 names(ffrefs) = cn
 ee = as(smlSet, "ExpressionSet")
 new("egtSet", assayData=assayData(ee),
	phenoData=phenoData(ee),
	experimentData=experimentData(ee),
	protocolData=protocolData(ee),
	featureData=featureData(ee),
	annotation=annotation(ee),
	gtlist=ffrefs,
	snpNames=sn)
}

setMethod("show", "egtSet", function(object) {
 callNextMethod()
 cat("Genotype information in a list of ff of length", length(object@gtlist), "\n")
 cat(" SNP data ncols:\n  ")
 cat(selectSome(sapply(object@gtlist, ncol)))
 cat("\n")
 cat(" SNP names for first element:\n  ")
 cat(selectSome(object@snpNames[[1]]))
 cat("\n")
})

setMethod("exprs", "egtSet", function(object) {
 assayData(object)$exprs
})

setMethod("smList", "egtSet", function(x) {
 x@gtlist
})

setGeneric("getSNP", function(egts, chrind, tok) standardGeneric("getSNP"))
setMethod("getSNP", c("egtSet", "ANY", "ScalarCharacter"), function(egts, chrind, tok) {
 sind = match(tok, egts@snpNames[[chrind]])
 smList(egts)[[chrind]][,sind,drop=FALSE]
})
