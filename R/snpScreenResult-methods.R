setMethod("show", "snpScreenResult", function(object) {
   cat("GGtools snpScreenResult for call:\n")
   print(object@call)
   if (!(object@fittertok  %in% c("fastAGM", "fastHET"))) {
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
  if (ssr@fittertok %in% c("fastAGM", "fastHET")) return(ssr[["pval"]])
  if (ssr@fittertok != "lm") stop("code is idiosyncratic for lm fits")
  ps = as.numeric(sapply(ssr, function(x) try(summary(x)$coef[2,4],silent=TRUE)))
}


plot_mlp = function (ssr, snpMeta, gchr = NULL, geneLocDF=NULL, ps = NULL, pch = 20, cex = 0.5, local = FALSE, plotf=smoothScatter) 
{
    if (ssr@fittertok %in% c("fastAGM", "fastHET"))
        ps = ssr[["pval"]]
    else if (ssr@fittertok != "lm") 
        stop("code is idiosyncratic for lm fits")
    if (is.null(ps)) 
        ps = as.numeric(sapply(ssr, function(x) try(summary(x)$coef[2, 
            4], silent = TRUE)))
    if (is.null(gchr)) {
        if (is.null(geneLocDF)) stop("if gchr omitted, please supply geneLocDF")
	if (!("chr" %in% names(geneLocDF))) stop("geneLocDF does not have chr column.")
    	gloc = geneLocDF[geneLocDF$gene == ssr@gene, ]
    	if (nrow(gloc) == 0) 
        	stop("can't determine gene's chromosome for title from geneLocDF, please supply as arg.")
    	gchr = gloc[1, "chr"]
    }
    bad = which(is.na(ps))
    if (!local) 
        XLIM = range(snpMeta[, "pos"])
    else {
        snpn = unlist(lapply(ssr, function(x) names(coef(x))[2]))
        XLIM = range(snpMeta[snpn, "pos"])
    }
    if (length(bad)>0) {
          ssr@locs = ssr@locs[-bad]
          ps = ps[-bad]
    }
    if (!is(snpMeta, "snpMetaWhole")) {
      plotf(ssr@locs, -log10(ps), xlab = paste("location on chromosome", 
          chromosome(snpMeta)), ylab = "-log10 p Ho:Bs=0", main = paste(ssr@gene, 
          "(chr", gchr, ")"), xlim = XLIM, pch = pch, cex = cex)
      if (!is.null(geneLocDF)) {
         for (i in 1:nrow(gloc)) {
             axis(3, at = gloc[i, "beg"], labels = FALSE, col = "green")
             axis(3, at = gloc[i, "end"], labels = FALSE, col = "red")
             }
      }
    }
    else if (is(snpMeta, "snpMetaWhole")) {
      plotf(ssr@locs, -log10(ps), xlab = paste("chromosome"),
          ylab = "-log10 p Ho:Bs=0", main = paste(ssr@gene, 
          "(chr", gchr, ")"), xlim = XLIM, pch = pch, cex = cex, axes=FALSE)
      axis(2)
      abline(v=snpMeta@chrbounds, col="gray")

      pts = c(0,snpMeta@chrbounds[-length(snpMeta@chrbounds)]) + diff(c(0,snpMeta@chrbounds))/2
      text ( pts, rep(-.1, length(pts) ) , snpMeta@chrlabs )

      if (!is.null(geneLocDF)) {
         for (i in 1:nrow(gloc)) {
             axis(3, at = gloc[i, "beg"], labels = FALSE, col = "green")
             axis(3, at = gloc[i, "end"], labels = FALSE, col = "red")
             }
      }
    }
    return(invisible(list(x = ssr@locs, y = -log10(ps))))
}
