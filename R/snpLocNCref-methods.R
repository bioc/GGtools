
setGeneric("getSnpData", function(pkg, ncdfRef) standardGeneric("getSnpData"))
setMethod("getSnpData", c("character", "character"), function(pkg, ncdfRef) {
  require(pkg, character.only=TRUE)
  locs = getSnpLocs(get(ncdfRef)(), gsub("_dbconn", "", ncdfRef))
  chr = getSnpChroms(get(ncdfRef)(), gsub("_dbconn", "", ncdfRef))
  list(chr=chr, cumloc=locs)
})
 
