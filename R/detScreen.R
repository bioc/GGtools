
detScreen = function(gge=c20GGceu, psn="206918_s_at", 
  chrmeta=chr20meta, chr="chr20", gran=50, gene="CPNE1", ...) {
opar = par()
cpn = regseq(gge, psn, 
  seq(1,ncol(gge@phenoData@pData),gran), chrmeta, ...)
par(mfrow=c(1,2))
plot(cpn$locs, -log10(cpn$pva), main=paste(psn,chr), xlab="position",
ylab="-log10 p Ho:B=0")
bot = which.min(cpn$pva)
ggrplot(gge, gene, names(bot))
par(opar)
invisible(list(bot=bot, cpn=cpn))
}

setClass("snpScreenResult", representation(call="call", 
    locs="numeric", chr="character", fittertok="character"), contains="list")

setGeneric("snpScreen", function(racExSet, snpMeta, gene, formTemplate, fitter, gran, ...)
    standardGeneric("snpScreen"))
setMethod("snpScreen", c("racExSet", "snpMeta", "genesym", "formula", "function", "numeric"),
    function(racExSet, snpMeta, gene, formTemplate=~., fitter=lm, gran=1, ...) {
      runTemplate = function(x,y) {
       z = as.character(x)
       z[2] = gsub("\\.", y, z[2])
       as.formula(z)
      }
      psid = getpsid( gene, annotation(racExSet) )
      y = exprs(racExSet)[psid,]
      outco = list(y)
      names(outco) = as.character(gene)
      nsnp = length(sn <- snpNames(racExSet))
      snpstodo = sn[ inuse <- seq(1, nsnp, gran) ]
      locs = snpMeta[sn, "pos", drop=TRUE]
      out = list()
      for (i in 1:length(snpstodo)) {
        if (options()$verbose == TRUE) {
           if (i %% 100 == 0) cat(i)
           }
        fm = runTemplate( formTemplate, snpstodo[i] ) 
        out[[i]] = try(oneFit(racExSet, outco, fm, fitter))
      }
      callsave = match.call()
      fittertok = deparse(substitute(fitter))
      new("snpScreenResult", call=callsave, 
        locs=locs, chr=chromosome(snpMeta), fittertok=fittertok, out)
})

setMethod("show", "snpScreenResult", function(object) {
   cat("GGtools snpScreenResult for call:\n")
   print(object@call)
   cat("There were", nf <- length(object), "attempted fits,\n")
   nerr = sum(sapply(object, function(x) inherits(x, "try-error")))
   cat("and", nf, "were successful.\n")
})

