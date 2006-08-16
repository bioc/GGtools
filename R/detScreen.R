
#detScreen = function(gge=c20GGceu, psn="206918_s_at", 
#  chrmeta=chr20meta, chr="chr20", gran=50, gene="CPNE1", ...) {
#opar = par()
#cpn = regseq(gge, psn, 
#  seq(1,ncol(gge@phenoData@pData),gran), chrmeta, ...)
#par(mfrow=c(1,2))
#plot(cpn$locs, -log10(cpn$pva), main=paste(psn,chr), xlab="position",
#ylab="-log10 p Ho:B=0")
#bot = which.min(cpn$pva)
#ggrplot(gge, gene, names(bot))
#par(opar)
#invisible(list(bot=bot, cpn=cpn))
#}

setClass("snpScreenResult", representation(call="call", gene="character", 
    locs="numeric", chr="character", fittertok="character"), contains="list")

setGeneric("snpScreen", function(racExSet, snpMeta, gene, formTemplate, fitter, gran, ...)
    standardGeneric("snpScreen"))
setMethod("snpScreen", c("racExSet", "snpMeta", "genesym", "formula", "function", "numeric"),
   function (racExSet, snpMeta, gene, formTemplate = ~., fitter = lm, 
      gran = 1, ...) 
  {
      runTemplate = function(x, y) {
          z = as.character(x)
          z[2] = gsub("\\.", y, z[2])
          as.formula(z)
      }
      psid = getpsid(gene, annotation(racExSet))
      y = exprs(racExSet)[psid, ]
      outco = list(y)
      names(outco) = as.character(gene)
      nsnp = length(sn <- snpNames(racExSet))
      snpstodo = sn[inuse <- seq(1, nsnp, gran)]
      if (any(is.na(snpstodo))) snpstodo = snpstodo[!is.na(snpstodo)]
      allpos = get("meta", snpMeta@meta)$pos
      allsn = rownames(get("meta", snpMeta@meta))
      names(allpos) = allsn
      snpstodo = intersect(snpstodo, allsn)
      locs = allpos[snpstodo]
      fittertok = deparse(substitute(fitter))
      callsave = match.call()
      out = list()
      if (fittertok == "fastAGM") {
        tmp = snps(racExSet)[snpstodo,]
        bad = apply(tmp,1,function(x) any(is.na(x)))
        if (any(bad)) {
           warning("some genotype results had missing values; associated SNPs are dropped completely in this version when fastAGM is used.")
           snpstodo = snpstodo[-which(bad)]
           }
        locs = allpos[snpstodo]
        ans = fastAGM(snps(racExSet)[snpstodo,], y)
        return(new("snpScreenResult", call=callsave, locs=locs, 
            chr=chromosome(snpMeta), fittertok=fittertok, gene=as.character(gene), ans))
      }
      for (i in 1:length(snpstodo)) {
          if (options()$verbose == TRUE) {
              if (i%%100 == 0) 
                  cat(i)
          }
          fm = runTemplate(formTemplate, snpstodo[i])
          out[[i]] = try(oneFit(racExSet, outco, fm, fitter))
      }
      names(out) = snpstodo
      new("snpScreenResult", call = callsave, locs = locs, chr = chromosome(snpMeta), 
        fittertok = fittertok, gene=as.character(gene), out)
})

setMethod("show", "snpScreenResult", function(object) {
   cat("GGtools snpScreenResult for call:\n")
   print(object@call)
   if (object@fittertok  != "fastAGM") {
     cat("There were", nf <- length(object), "attempted fits,\n")
     nerr = sum(sapply(object, function(x) inherits(x, "try-error")))
     cat("and", nf-nerr, "were successful.\n")
   }
   else {
     cat("There were", nf <- length(object[[4]]), "attempted fits,\n")
     nerr = sum(is.na(object[[4]])) + sum(is.nan(object[[4]]))
     cat("and", nf-nerr, "were successful.\n")
   }
})

extract_p = function(ssr) {
  if (ssr@fittertok == "fastAGM") return(ssr[["pval"]])
  if (ssr@fittertok != "lm") stop("code is idiosyncratic for lm fits")
  ps = as.numeric(sapply(ssr, function(x) try(summary(x)$coef[2,4],silent=TRUE)))
}

plot_mlp <- function (ssr, snpMeta, ps=NULL, pch=20, cex=1) 
{
    if (ssr@fittertok == "fastAGM") ps = ssr[["pval"]]
    else if (ssr@fittertok != "lm") 
        stop("code is idiosyncratic for lm fits")
    if (is.null(ps)) ps = as.numeric(sapply(ssr, function(x) try(summary(x)$coef[2, 
        4], silent = TRUE)))
    if (length(ssr@locs) > 200) plotf = .smoothScatter2
     else plotf = plot
    bad = which(is.na(ps))
    plotf(ssr@locs[-bad], -log10(ps[-bad]), xlab = "chromosomal location", 
        ylab = "-log10 p Ho:Bs=0", main = ssr@gene, xlim=range(snpMeta[,"pos"]),        pch=pch, cex=cex)
    data(geneLocs)
    x = geneLocs[geneLocs$gene == ssr@gene, ]
    if (nrow(x) == 0) 
        return(invisible(NULL))
    for (i in 1:nrow(x)) {
        axis(3, at = x[i, "beg"], labels = FALSE, col = "green")
        axis(3, at = x[i, "end"], labels = FALSE, col = "red")
    }
    return(invisible(NULL))
}

