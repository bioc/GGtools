setClassUnion("integerOrNULL", c("integer", "NULL"))

#function (smpack, rhs = ~1, nperm = 2, folderstem = "cisScratch", 
#    radius = 50000, shortfac = 100, chrnames = "22", smchrpref = "", 
#    gchrpref = "", schrpref = "ch", geneApply = lapply, geneannopk = "illuminaHumanv1.db", 
#    snpannopk = snplocsDefault(), smFilter = function(x) nsFilter(MAFfilter(x, 
#        lower = 0.05), var.cutoff = 0.9), exFilter = function(x) x, 
#    keepMapCache = FALSE, SSgen = GGBase::getSS, excludeRadius = NULL, 
#    estimates = FALSE, ...) 

setClass("CisConfig", representation(
  smpack = "character",
  rhs = "formula",
  nperm = "integer",
  folderStem = "character",
  radius = "integer",
  MAFlb = "numeric",
  shortfac = "integer",
  chrnames = "character",
  smchrpref = "character",
  gchrpref = "character",
  schrpref = "character",
  geneApply = "function",
  geneannopk = "character",
  snpannopk = "character",
  smFilter = "function",
  exFilter = "function",
  keepMapCache = "logical",
  SSgen = "function",
  excludeRadius = "integerOrNULL",
  estimates = "logical"))

setMethod("show", "CisConfig", function(object) {
 cat("CisConfig instance.  Key parameters:\n")
 cat("smpack = ", smpack(object), "; chrnames = ", chrnames(object), "\n")
 cat("nperm = ", nperm(object), "; radius = ", radius(object), "\n====\n")
 cat("Configure using \n")
 print(paste0(slotNames(new("CisConfig")), "<-"))
})

setMethod("initialize", "CisConfig", function(.Object) {
 .Object@smpack = "GGdata"
 .Object@rhs = ~1
 .Object@nperm = 2L
 .Object@folderStem = "cisScratch"
  .Object@radius = 50000L
  .Object@MAFlb = .025
  .Object@shortfac = 100L
  .Object@chrnames = "22"
  .Object@smchrpref = ""
  .Object@gchrpref = ""
  .Object@schrpref = "ch"
  .Object@geneApply = lapply
  .Object@geneannopk = "illuminaHumanv1.db"
  .Object@snpannopk = snplocsDefault()
  .Object@smFilter = function(x) nsFilter(MAFfilter(x, lower=.05), var.cutoff=.8)
  .Object@exFilter = force
  .Object@keepMapCache = FALSE
  .Object@SSgen = GGBase::getSS
  .Object@excludeRadius = 0L
  .Object@estimates = TRUE
  .Object
})

setGeneric("smpack", function(x) standardGeneric("smpack"))
setMethod("smpack", "CisConfig", function(x) x@smpack)
setGeneric("smpack<-", function(object, value) standardGeneric("smpack<-"))
setMethod("smpack<-", c("CisConfig", "character"), function(object, value) {object@smpack <- value; object})
setGeneric("smFilter", function(x) standardGeneric("smFilter"))
setMethod("smFilter", "CisConfig", function(x) x@smFilter)
setGeneric("smFilter<-", function(object, value) standardGeneric("smFilter<-"))
setMethod("smFilter<-", c("CisConfig", "function"), function(object, value) {object@smFilter <- value; object})
setGeneric("rhs", function(x) standardGeneric("rhs"))
setMethod("rhs", "CisConfig", function(x) x@rhs)
setGeneric("rhs<-", function(object, value) standardGeneric("rhs<-"))
setMethod("rhs<-", c("CisConfig", "function"), function(object, value) {object@rhs <- value; object})
setGeneric("nperm", function(x) standardGeneric("nperm"))
setMethod("nperm", "CisConfig", function(x) x@nperm)
setGeneric("nperm<-", function(object, value) standardGeneric("nperm<-"))
setMethod("nperm<-", c("CisConfig", "integer"), function(object, value) {object@nperm <- value; object})
setGeneric("folderStem", function(x) standardGeneric("folderStem"))
setMethod("folderStem", "CisConfig", function(x) x@folderStem)
setGeneric("folderStem<-", function(object, value) standardGeneric("folderStem<-"))
setMethod("folderStem<-", c("CisConfig", "character"), function(object, value) {object@folderStem <- value; object})
setGeneric("radius", function(x) standardGeneric("radius"))
setMethod("radius", "CisConfig", function(x) x@radius)
setGeneric("radius<-", function(object, value) standardGeneric("radius<-"))
setMethod("radius<-", c("CisConfig", "integer"), function(object, value) {object@radius <- value; object})
setGeneric("MAFlb", function(object, value) standardGeneric("MAFlb"))
setMethod("MAFlb", c("CisConfig"), function(object) {object@MAFlb})
setGeneric("MAFlb<-", function(object, value) standardGeneric("MAFlb<-"))
setMethod("MAFlb<-", c("CisConfig", "integer"), function(object, value) {object@MAFlb <- value; object})
setGeneric("shortfac", function(x) standardGeneric("shortfac"))
setMethod("shortfac", "CisConfig", function(x) x@shortfac)
setGeneric("shortfac<-", function(object, value) standardGeneric("shortfac<-"))
setMethod("shortfac<-", c("CisConfig", "integer"), function(object, value) {object@shortfac <- value; object})
setGeneric("chrnames", function(x) standardGeneric("chrnames"))
setMethod("chrnames", "CisConfig", function(x) x@chrnames)
setGeneric("chrnames<-", function(object, value) standardGeneric("chrnames<-"))
setMethod("chrnames<-", c("CisConfig", "character"), function(object, value) {object@chrnames <- value; object})
setGeneric("smchrpref", function(x) standardGeneric("smchrpref"))
setMethod("smchrpref", "CisConfig", function(x) x@smchrpref)
setGeneric("smchrpref<-", function(object, value) standardGeneric("smchrpref<-"))
setMethod("smchrpref<-", c("CisConfig", "character"), function(object, value) {object@smchrpref <- value; object})
setGeneric("gchrpref", function(x) standardGeneric("gchrpref"))
setMethod("gchrpref", "CisConfig", function(x) x@gchrpref)
setGeneric("gchrpref<-", function(object, value) standardGeneric("gchrpref<-"))
setMethod("gchrpref<-", c("CisConfig", "character"), function(object, value) {object@gchrpref <- value; object})
setGeneric("schrpref", function(x) standardGeneric("schrpref"))
setMethod("schrpref", "CisConfig", function(x) x@schrpref)
setGeneric("schrpref<-", function(object, value) standardGeneric("schrpref<-"))
setMethod("schrpref<-", c("CisConfig", "character"), function(object, value) {object@schrpref <- value; object})
setGeneric("geneApply", function(x) standardGeneric("geneApply"))
setMethod("geneApply", "CisConfig", function(x) x@geneApply)
setGeneric("geneApply<-", function(object, value) standardGeneric("geneApply<-"))
setMethod("geneApply<-", c("CisConfig", "function"), function(object, value) {object@geneApply <- value; object})
setGeneric("geneannopk", function(x) standardGeneric("geneannopk"))
setMethod("geneannopk", "CisConfig", function(x) x@geneannopk)
setGeneric("geneannopk<-", function(object, value) standardGeneric("geneannopk<-"))
setMethod("geneannopk<-", c("CisConfig", "character"), function(object, value) {object@geneannopk <- value; object})
setGeneric("snpannopk", function(x) standardGeneric("snpannopk"))
setMethod("snpannopk", "CisConfig", function(x) x@snpannopk)
setGeneric("snpannopk<-", function(object, value) standardGeneric("snpannopk<-"))
setMethod("snpannopk<-", c("CisConfig", "character"), function(object, value) {object@snpannopk <- value; object})
setGeneric("exFilter", function(x) standardGeneric("exFilter"))
setMethod("exFilter", "CisConfig", function(x) x@exFilter)
setGeneric("exFilter<-", function(object, value) standardGeneric("exFilter<-"))
setMethod("exFilter<-", c("CisConfig", "function"), function(object, value) {object@exFilter <- value; object})
setGeneric("keepMapCache", function(x) standardGeneric("keepMapCache"))
setMethod("keepMapCache", "CisConfig", function(x) x@keepMapCache)
setGeneric("keepMapCache<-", function(object, value) standardGeneric("keepMapCache<-"))
setMethod("keepMapCache<-", c("CisConfig", "logical"), function(object, value) {object@keepMapCache <- value; object})
setGeneric("SSgen", function(x) standardGeneric("SSgen"))
setMethod("SSgen", "CisConfig", function(x) x@SSgen)
setGeneric("SSgen<-", function(object, value) standardGeneric("SSgen<-"))
setMethod("SSgen<-", c("CisConfig", "function"), function(object, value) {object@SSgen <- value; object})
setGeneric("excludeRadius", function(x) standardGeneric("excludeRadius"))
setMethod("excludeRadius", "CisConfig", function(x) x@excludeRadius)
setGeneric("excludeRadius<-", function(object, value) standardGeneric("excludeRadius<-"))
setMethod("excludeRadius<-", c("CisConfig", "integerOrNULL"), function(object, value) {object@excludeRadius <- value; object})
setGeneric("estimates", function(x) standardGeneric("estimates"))
setMethod("estimates", "CisConfig", function(x) x@estimates)
setGeneric("estimates<-", function(object, value) standardGeneric("estimates<-"))
setMethod("estimates<-", c("CisConfig", "logical"), function(object, value) {object@estimates <- value; object})
