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
   snplocs="GRanges", generanges="GRanges", radiusUsed="numeric"))

setClass("cwBestCis", contains="RangedData")
setClass("mcwBestCis", representation(scoregr = "GRanges",
 allperm="numeric", extra="ANY", chromUsed="ANY", theCall="call"))


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

