
setGeneric("getSnpData", function(x) standardGeneric("getSnpData"))
setMethod("getSnpData", "character", function(x) {
  oo = open.ncdf(x)
  on.exit(close.ncdf(oo))
  list(chr=get.var.ncdf(oo, "chr"), cumloc=get.var.ncdf(oo, "cumloc"))
})
 
