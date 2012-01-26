#gsetFmla2FmlaList = function(fm) {
##
## this function takes a formula with a GeneSet on the lhs
## and returns a list of gene-specific formulas
##
#.Deprecated("eqtlTests", , msg="please use eqtlTests for all GGtools analyses")
# flist = as.list(fm)
# gs = eval(flist[[2]])
# if (length(flist[[3]]) > 1) pred = paste(flist[[3]][-1], collapse="+")
# else if (length(flist[[3]]) == 1) pred = flist[[3]]
# if (!is(gs, "GeneSet")) stop("needs GeneSet instance in response position")
# wrapg = function(x) paste("genesym(\"", x, "\")", sep="")
# wrapex = function(x) paste("probeId(\"", x, "\")", sep="")
# toks = GSEABase::geneIds(gs)
# idty = geneIdType(gs)
# if (is(idty, "SymbolIdentifier"))
#   resps = wrapg(toks)
# else if (is(idty, "AnnotationIdentifier"))
#   resps = wrapex(toks)
# else stop("geneIdType for gene set must be either SymbolIdentifier or AnnotationIdentifier")
# lapply(paste(resps, pred, sep="~"), formula)
#}
#
   
