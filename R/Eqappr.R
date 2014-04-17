
.validEQA = function(object) isTRUE(all(c("sens", "unpruned", "pruned") %in% names(object)))

comlist = function(x) paste(x, collapse=", ")

setClass("ApprMeta", representation(pop="character", hetcorrect="character"))

setClass("ApprSens", representation(meta="ApprMeta", sens="list"))

setClass("ApprRes", representation(meta="ApprMeta", 
    prunedtype="character", coeflist="list",
    tabs="list", outs="list", dims="list"), contains="list")

setClass("EqAppr", representation(meta="ApprMeta", sens="ApprSens",
    pruned="ApprRes", unpruned="ApprRes"))

setMethod("show", "EqAppr", function(object) {
 cat("eQTL appraisal container, population ", sQuote(object@meta@pop),
       " (het. correct.: ", sQuote(object@meta@hetcorrect),  ")\n", sep="")
 cat("  sensitivities available for each of", comlist(names(object@sens@sens)), "\n")
 cat("  appraisal components:\n")
 show(object@unpruned) 
 show(object@pruned) 
 cat(" For each ApprRes model component, elements \n")
  cat("   ", comlist(names(object@pruned@outs[[1]])), "\n")
 cat(" are available.\n")
})

#setGeneric("Ntest", function(x) standardGeneric("Ntest"))
#setGeneric("Ntrain", function(x) standardGeneric("Ntrain"))
#setMethod("Ntest", "Eqappr", function(x) {
#   ans = sapply(x[-1], function(x) x$dimtest[1])
#   names(ans) = names(x)[-1]
#   ans
#})
#
#setMethod("Ntrain", "Eqappr", function(x) {
#   ans = sapply(x[-1], function(x) x$dimtrain[1])
#   names(ans) = names(x)[-1]
#   ans
#})

#setGeneric("unpruned", function(x) standardGeneric("unpruned"))
#setMethod("unpruned", "Eqappr", function(x) x$unpruned)
#setGeneric("pruned", function(x) standardGeneric("pruned"))
#setMethod("pruned", "Eqappr", function(x) x$pruned)
#

buildApprRes = function( inob, inmeta, type="pruned" ) {
 new("ApprRes", meta=inmeta, prunedtype=type, coeflist=inob[[type]]$coeflist,
   tabs=inob[[type]]$tabs,
   outs=inob[[type]]$outs,
   dims=list(dtab=inob[[type]]$dimdtab,
   dimtest=inob[[type]]$dimtest,
   dimtrain=inob[[type]]$dimtrain)) }

#> names(ceu_pc10[[2]])
#[1] "coeflist" "tabs"     "outs"     "dimdtab"  "dimtest"  "dimtrain"

setMethod("show", "ApprRes", function(object) {
 cat(" ApprRes instance of type", sQuote(object@prunedtype), "\n")
 cat("  ", length(object@coeflist), "models checked, including\n")
 cat("    ", selectSome(names(object@coeflist)), "\n")
 })

allBinnedPredictions = function(apres) sapply(apres@tabs, force)
allCalibProportions = function(apres) {
  ans = sapply(apres@coeflist, function(x)x[,1])
  rownames(ans) = sub(".*\\(", "(", rownames(ans))
  ans }
allTruePos = function(apres) {
  ans = allCalibProportions(apres)*allBinnedPredictions(apres)
  rownames(ans) = sub(".*\\(", "(", rownames(ans))
  ans }
getAUCs = function(apres) unlist(sapply(apres@outs, function(x)x$auc@y.values))

#library(biglm)
#library(ROCR)

setGeneric("getSens", function(x, ...) standardGeneric("getSens"))
setMethod("getSens", "EqAppr", function(x, ...) x@sens)
setGeneric("getPruned", function(x, ...) standardGeneric("getPruned"))
setMethod("getPruned", "EqAppr", function(x, ...) x@pruned)
setGeneric("getUnpruned", function(x, ...) standardGeneric("getUnpruned"))
setMethod("getUnpruned", "EqAppr", function(x, ...) x@unpruned)
setGeneric("getNrec", function(x, ...) standardGeneric("getNrec"))
setMethod("getNrec", "ApprRes", function(x, ...) sapply(x@dims, "[", 1))
setGeneric("getModnames", function(x)standardGeneric("getModnames"))
setMethod("getModnames", "EqAppr", function(x) names(x@unpruned@coeflist))

buildEqAppr = function(inobj, pop, hetcorrect) {
 meta = new("ApprMeta", pop=pop, hetcorrect=hetcorrect)
 sens = new("ApprSens", meta=meta, sens=inobj$sens )
 pruned = buildApprRes( inobj, meta, type="pruned" )
 unpruned = buildApprRes( inobj, meta, type="unpruned" )
 new("EqAppr", meta=meta, pruned=pruned, unpruned=unpruned,
   sens=sens)
}

#load("ceu_pc15.rda")
#aresC15 = buildEqAppr( ceu_pc15, "CEU", "15 PC")
#load("ceu_pc5.rda")
#aresC5 = buildEqAppr( ceu_pc5, "CEU", "5 PC")
#load("ceu_pc10.rda")
#aresC10 = buildEqAppr( ceu_pc10, "CEU", "10 PC")
#
#load("yri_pc15.rda")
#aresY15 = buildEqAppr( yri_pc15, "YRI", "15 PC")
#load("yri_pc5.rda")
#aresY5 = buildEqAppr( yri_pc5, "YRI", "5 PC")
#load("yri_pc10.rda")
#aresY10 = buildEqAppr( yri_pc10, "YRI", "10 PC")

maxAtFDR = function( ea, fdtag = "at_0.1", type="probes" ) {
 max(getSens( ea )@sens[[type]][, fdtag])
}

optAUC = function( ea, type="pruned" ) {
 ap = getUnpruned(ea)
 if (type=="pruned") ap = getPruned(ea)
 max(getAUCs(ap))
}

.calfig = function (colist, tabs, ind = 10, hfudgetxt = 0.0155, tickend = 0.16, 
    tickgap = 0.02, ylimin = c(-0.01, 0.16), xlimin = c(-0.01, 
        0.16), fraccex = 0.8, fuselast = 0, add=FALSE) 
{
    getmidpts = function(coefrn) {
        lima = gsub(".*\\(", "", coefrn)
        limb = gsub("]", "", lima)
        limc = strsplit(limb, ",")
        limd = sapply(limc, as.numeric)
        apply(limd, 2, mean)
    }
    midcuts = getmidpts(rownames(colist[[ind]]))
  if (!add) {
    plot(0, 0, ylim = ylimin, xlim = xlimin, pch = " ", xlab = "predicted", 
        ylab = "empirical", axes = FALSE)
    axis(1, at = seq(0, tickend, tickgap))
    axis(2, at = seq(0, tickend, tickgap))
    }
    fracs = paste(nums <- round(tabs[[ind]] * colist[[ind]][, 
        1]), tabs[[ind]], sep = "/")
    if (fuselast == 0) {
        points(midcuts, colist[[ind]][, 1], col=ifelse(add, "black", "red"))
        if (!add) text(midcuts + hfudgetxt, colist[[ind]][, 1], labels = fracs, 
            cex = fraccex)
    }
    else {
        l.m = length(midcuts)
        idrop = (l.m - fuselast + 1):l.m
        midcuts = c(midcuts[-idrop], mean(midcuts[idrop]))
        ndrp = tabs[[ind]][idrop]
        cdrp = colist[[ind]][idrop, 1]
        hdrp = cdrp * ndrp
        newfrac = sum(hdrp)/sum(ndrp)
        newfracc = paste(round(sum(hdrp), 0), sum(ndrp), sep = "/")
        fracs = c(fracs[-idrop], newfracc)
        newco = colist[[ind]][-idrop, 1]
        newco = c(newco, newfrac)
        points(midcuts, newco, col=ifelse(add, "black", "red"))
        if (!add) text(midcuts + hfudgetxt, newco, labels = fracs, cex = fraccex)
    }
    abline(0, 1, lty = 2, col = "gray")
}

setGeneric("calfig", function(x, ind, ...) standardGeneric("calfig"))
setMethod("calfig", c("EqAppr", "character"), function(x, ind, ...) {
  tmp = getUnpruned(x) 
  colist = tmp@coeflist
  tablist = tmp@tabs
 .calfig(colist, tablist, ind, ...)
})
