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
  genome = "character",
  excludeRadius = "integerOrNULL",
  estimates = "logical",
  extraProps="function", useME="logical", MEpvot="numeric"))

setMethod("show", "CisConfig", function(object) {
 cat("CisConfig instance; genome ", genome(object),".  Key parameters:\n")
 cat("smpack = ", smpack(object), "; chrnames = ", chrnames(object), "\n")
 cat("nperm = ", nperm(object), "; radius = ", radius(object), "\n====\n")
 cat("Configure using \n")
 print(paste0(slotNames(new("CisConfig")), "<-"))
})

setMethod("initialize", "CisConfig", function(.Object) {
 .Object@smpack = "GGdata"
 .Object@rhs = ~1
 .Object@genome = "hg19"
 .Object@nperm = 3L
 .Object@folderStem = "cisScratch"
  .Object@radius = 50000L
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
  .Object@keepMapCache = TRUE
  .Object@SSgen = GGBase::getSS
  .Object@excludeRadius = 0L
  .Object@estimates = TRUE
  .Object@extraProps = force
  .Object@useME = FALSE
  .Object@MEpvot = .5
  .Object
})

setGeneric("smpack", function(x) standardGeneric("smpack"))
setMethod("genome", "CisConfig", function(x) x@genome)
setMethod("genome<-", "CisConfig", function(x, value) {x@genome = value; x})
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
setMethod("rhs<-", c("CisConfig", "formula"), function(object, value) {object@rhs <- value; object})
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
setGeneric("extraProps", function(x) standardGeneric("extraProps"))
setGeneric("extraProps<-", function(object,value) standardGeneric("extraProps<-"))
setMethod("extraProps", "CisConfig", function(x) x@extraProps)
setMethod("extraProps<-", c("CisConfig", "function"), 
   function(object, value) {object@extraProps <- value; object})


setClass("TransConfig", contains="CisConfig",
 representation(snpchr="character", gbufsize="integer",
  batchsize="integer"))

setGeneric("snpchr", function(x) standardGeneric("snpchr"))
setGeneric("snpchr<-", function(object,value) standardGeneric("snpchr<-"))
setMethod("snpchr", "TransConfig", function(x) x@snpchr)
setMethod("snpchr<-", c("TransConfig", "character"), function(object,value) {object@snpchr <- value; object})

setGeneric("gbufsize", function(x) standardGeneric("gbufsize"))
setGeneric("gbufsize<-", function(object,value) standardGeneric("gbufsize<-"))
setMethod("gbufsize", "TransConfig", function(x) x@gbufsize)
setMethod("gbufsize<-", c("TransConfig", "integer"), 
   function(object, value) {object@gbufsize <- value; object})

setGeneric("batchsize", function(x) standardGeneric("batchsize"))
setGeneric("batchsize<-", function(object,value) standardGeneric("batchsize<-"))
setMethod("batchsize", "TransConfig", function(x) x@batchsize)
setMethod("batchsize<-", c("TransConfig", "integer"), 
   function(object, value) {object@batchsize <- value; object})

setMethod("initialize", "TransConfig", function(.Object) {
 .Object = callNextMethod()
 .Object@radius = 100000L
 .Object@chrnames = as.character(1:22)
 .Object@gbufsize = 20L
 .Object@batchsize = 200L
 .Object
})

setMethod("show", "TransConfig", function (object) 
{
    cat("TransConfig instance.  Key parameters:\n")
    cat("smpack = ", smpack(object), "; snpchr = ", snpchr(object), "; chrnames = ", selectSome(chrnames(object)), 
        "\n")
    cat("nperm = ", nperm(object), "; radius = ", radius(object), 
        "\n====\n")
    cat("Configure using \n")
    print(paste0(slotNames(new("TransConfig")), "<-"))
})

ivector <- function(x, ...) {
 i <- 1
 it <- idiv(length(x), ...)

 nextEl <- function() {
   n <- nextElem(it)
   ix <- seq(i, length=n)
   i <<- i + n
   x[ix]
 }

 obj <- list(nextElem=nextEl)
 class(obj) <- c('ivector', 'abstractiter', 'iter')
 obj
}

add878 = function(ans) {
  data(hmm878)
  ac = as.character
  eqr = GRanges(ac(seqnames(ans)), IRanges(ans$snplocs, width=1))
  fo = findOverlaps(eqr, hmm878)
  chromcat878 = factor(rep("none", length(ans)), levels=c(unique(hmm878$name), "none"))
  chromcat878[ queryHits(fo) ] = factor(hmm878$name[subjectHits(fo)])
  ans$chromcat878 = chromcat878
  ans
}

inflammFilter = function(gwtagger) {
  require(gwascat)
# gwrngs in scope
  allt = gwrngs$Disease.Trait
  infinds = grep("rheumatoid|inflamm|crohn|lupus|multiple sclero|type 1 diabetes",
     allt, ignore.case=TRUE)
  gwtagger[ which(overlapsAny( gwtagger, gwrngs[infinds]) | gwtagger$baseid %in% gwrngs[infinds]$SNPs) ]
}

addgwhit = function(ans, traitFilter=force, vname="isgwashit") {
    if (require(gwascat)) {
    data(gwastagger)
    ac = as.character
    if (is(ans, "data.table")) seqn = as.character(ans$seqnames)
    else if (inherits(ans, "GRanges")) seqn = ac(seqnames(ans))
    else stop("ans not data.table or GRanges derivative")
    eqr = GRanges(seqn, IRanges(ans$snplocs, width=1))
    gwt = traitFilter(gwastagger)
    isgwashit = 1*(overlapsAny(eqr, gwt) | ans$snp %in% gwt$tagid) # allow match by loc or name
    if (is(ans, "data.table")) 
       ans[[vname]] = isgwashit
    else mcols(ans)[,vname] = isgwashit
    }
  else warning("gwascat not available; returning ans unaltered")
  ans
}

get_probechunks = function(smpack="yri1kgv", chrpref="chr", chunksize=250,
   allc=21:22) {
  sm = getSS(smpack, paste0(chrpref, "22")) # illustrative, source of probeids
  allp = featureNames(sm)
  ganno = annotation(sm)
  require(ganno, character.only=TRUE)
  cmap = select(get(ganno), keytype="PROBEID", keys=allp, columns="CHR")
  cmap = cmap[which(cmap$CHR %in% as.character(allc)),]
  byc = split(cmap$PROBEID, cmap$CHR)
  alli = lapply(byc, ivector, chunkSize=chunksize) # get nice balance
  lapply(alli, as.list) # materialize
}

buildConfList = function( baseconf, chunksize = 100, chromToDo=1:22 ) {
  smpack = smpack(baseconf)
  pchunks = get_probechunks( smpack=smpack, 
      chrpref=smchrpref(baseconf), chunksize=chunksize,
      allc = chromToDo )
  nel = sum(clen <- sapply(pchunks, length))
  cnames = rep(names(pchunks), clen)
  configList = vector("list", nel)
  plist = unlist(pchunks, recursive=FALSE)
  for (i in 1:nel) {
    tmp = baseconf
    z = function() function(x) smFilter(baseconf)(x)[probeId(
      intersect(featureNames(x),pl)),]
    smFilter(tmp) = z()  # must skirt lazy evaluation
    environment(smFilter(tmp))$pl = plist[[i]]
    chrnames(tmp) = as.character(cnames[i])
    folderStem(tmp) = paste0(folderStem(tmp), "_", cnames[i], "_",
        plist[[i]][1])
    configList[[i]] = tmp
    }
  configList
}
