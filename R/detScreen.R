
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

setGeneric("snpScreen", function(racExSet, snpMeta, gene, formTemplate, fitter, gran, ...)
    standardGeneric("snpScreen"))
setMethod("snpScreen", c("racExSet", "snpMeta", "genesym", "formula", "function", "numeric"),
    function(racExSet, snpMeta, gene, formTemplate=~., fitter=lm, gran=1, ...) {
      runTemplate = function(x,y) { 
       # assuming that a template of the form ~. or ~factor(.) is supplied
       # this replaces . by string y
          uu = as.name(gsub("\\.", x[[2]], y))
          x[[2]] = substitute(uu)
          x
      }
      psid = getpsid( gene, annotation(racExSet) )
      y = exprs(racExSet)[psid,]
      outco = list(y)
      names(outco) = as.character(gene)
      nsnp = length(sn <- snpNames(racExSet))
      snpstodo = sn[ seq(1, nsnp, gran) ]
      out = list()
      for (i in 1:length(snpstodo)) {
        fm = runTemplate( formTemplate, snpstodo[i] ) 
        out[[i]] = oneFit(racExSet, outco, fm, fitter)
      }
      out
})

