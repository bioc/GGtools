
setGeneric("getSnpData", function(pkg, ncdfRef) standardGeneric("getSnpData"))
setMethod("getSnpData", c("character", "character"), function(pkg, ncdfRef) {
  require(pkg, character.only=TRUE)
  x = get( ncdfRef, paste("package:", pkg, sep=""))
  list(chr=get.var.ncdf(x, "chr"), cumloc=get.var.ncdf(x, "cumloc"))
})
 
