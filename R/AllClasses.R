setOldClass("sessionInfo")
setClass("SessionInfo", contains="sessionInfo")

setClass("gwSnpScreenResult", contains="list",
   representation(gene="character", psid="character", annotation="character", chrnum="ANY",
      testType="character", call="call",
        sessionInfo="SessionInfo")) #, modFmla="formula"))

setOldClass("ff_matrix")
setOldClass("ff_array")

chkeman = function(object){
# eqtlTestsManager validity test
 return(TRUE)
}

chkeeman = chkeman  # eventually enlarge for estimates object


setClass("eqtlTestsManager",
 representation(fffile="ff_matrix", call="call", sess="ANY",
        exdate="ANY", shortfac="numeric", geneanno="character", df="numeric",
        summaryList="list"),
        validity=chkeman)

setClass("eqtlEstimatesManager", contains="eqtlTestsManager",
        validity=chkeeman)

setClass("cisMap", representation(namelist="list",
   snplocs="GRanges", generanges="GRanges", radiusUsed="numeric", excludeRadius="numeric"))

setClass("cwBestCis", contains="GRanges")
setClass("mcwBestCis", representation(scoregr = "GRanges",
 allperm="numeric", extra="ANY", chromUsed="ANY", theCall="call", smFilter="function", nperm="numeric", globalMap="environment",
 testCount="numeric"))
setClass("allSigCis", representation(fulllist = "GRanges", bestcis="mcwBestCis",
 chromUsed="ANY", theCall="call"))
#
# could compute an approximate FDR for all elements of the allBestCis using the
# allperm component of the mcwBestCis component
#


setClass("transManager", representation(base="list"))

setMethod("show", "transManager", function(object){
 basel = object@base
 cat("transManager instance, created", basel$date, "\n", sep=" ")
 cat("dimension of scores component:\n")
 cat(" number of loci checked: ", nrow(basel$scores),
   "; genes retained: ", ncol(basel$scores), "\n", sep="")
 cat("the call was:\n")
 print(basel$call)
})

.combine2 = function( mcw1, mcw2 ) {
# rudimentary combination
 thecall = match.call()
 nperm = mcw1@nperm  # FIXME NEED TO PULL FROM A NEW SLOT!
 if (nperm != mcw2@nperm) stop("two inputs have different nperm fields")
 obs = suppressWarnings(c(mcw1@scoregr, mcw2@scoregr)) # uses different seq sets
 alls = c( mcw1@allperm, mcw2@allperm )
 fdrs = sapply(elementMetadata(obs)$score, function(x) (sum(alls>=x)/nperm)/sum(elementMetadata(obs)$score>=x))
 elementMetadata(obs)$fdr = fdrs
 obs = obs[ order(elementMetadata(obs)$fdr) , ]
 new("mcwBestCis", scoregr = obs, allperm = alls, theCall = thecall, chromUsed = c(mcw1@chromUsed, mcw2@chromUsed),
   extra = list(mcw1@smFilter, mcw2@smFilter), nperm=nperm )
}

combine2 = function (mcw1, mcw2, doFDR=TRUE) {
#
# unexported
#
    thecall = match.call()
    nperm = mcw1@nperm
    if (nperm != mcw2@nperm) 
        stop("two inputs have different nperm fields")
    obs = suppressWarnings(c(mcw1@scoregr, mcw2@scoregr))
    alls = c(mcw1@allperm, mcw2@allperm)
if (doFDR) {
    fdrs = sapply(elementMetadata(obs)$score, function(x) (sum(alls >= 
        x)/nperm)/sum(elementMetadata(obs)$score >= x))
    elementMetadata(obs)$fdr = fdrs
    obs = obs[order(elementMetadata(obs)$fdr), ]
}
    new("mcwBestCis", scoregr = obs, allperm = alls, theCall = thecall, 
        chromUsed = c(mcw1@chromUsed, mcw2@chromUsed), extra = list(mcw1@smFilter, 
            mcw2@smFilter), nperm = nperm)
}

setClass("mcwAllCis", representation(obs="GRanges", perms="list",
    theCall="call"))
setMethod("show", "mcwAllCis", function(object) {
 cat("mcwAllCis instance.\n")
 show(object@obs)
})
