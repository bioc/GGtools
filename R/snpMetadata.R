
setGeneric("getSnpMetaConn", function(x) standardGeneric("getSnpMetaConn"))

setMethod("getSnpMetaConn", "smlSet", function(x) {
  get(slot(x, "snpLocRef"))()
})

setMethod("getSnpMetaConn", "cwSnpScreenResult", function(x) {
  get(slot(x, "snpLocExtRef"))()
})

setGeneric("getSnpMetadata", function(x, sel)
  standardGeneric("getSnpMetadata"))

setMethod("getSnpMetadata", c("smlSet", "missing"),
  function(x, sel) {
  con = getSnpMetaConn(x)
  tab = dbListTables(con)
  que = paste("select * from ", tab, sep="")
  dbGetQuery(con, que)
})

setMethod("getSnpMetadata", c("smlSet", "chrnum"),
  function(x, sel) {
  con = getSnpMetaConn(x)
  tab = dbListTables(con)
  que = paste("select * from ", tab, " where chrnum = ", 
    as.character(sel), sep="")
  dbGetQuery(con, que)
})

setMethod("getSnpMetadata", c("smlSet", "rsNum"),
  function(x, sel) {
  rsids = gsub("rs", "", as.character(sel))
  if (length(rsids) < 1) stop("no valid RS numbers supplied")
  if (length(rsids)>1) incl = paste("(", paste(rsids, collapse=","), ")", sep="")
  else incl = paste("(", rsids, ")", sep="")
  con = getSnpMetaConn(x)
  tab = dbListTables(con)
  que = paste("select * from ", tab, " where rsid in ", 
    incl, sep=" ")
  dbGetQuery(con, que)
})

setGeneric("getSnpOffsets", function(x) standardGeneric("getSnpOffsets"))

setMethod("getSnpOffsets", "smlSet", function(x) {
  con = getSnpMetaConn(x)
  tab = dbListTables(con)
  que = paste("select * from ", tab, " limit 24", sep="")
  dbGetQuery(con, que)[,"offsets"]
})
 
setMethod("getSnpOffsets", "cwSnpScreenResult", function(x) {
  con = getSnpMetaConn(x)
  tab = dbListTables(con)
  que = paste("select * from ", tab, " limit 24", sep="")
  dbGetQuery(con, que)[,"offsets"]
})
