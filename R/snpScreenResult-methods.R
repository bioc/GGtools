setMethod("annotation", "snpScreenResult", function(object) object@annotation)

setMethod("show", "snpScreenResult", function(object) {
   cat("GGtools snpScreenResult for call:\n")
   print(object@call)
     cat("There were", nf <- length(object[[4]]), "attempted fits,\n")
     cat("and", length(na.omit(object[[4]])) , "were successful.\n")
})

extract_p = function(ssr) {
  if (is(ssr@fitter, "GGfitter")) return(ssr[["pval"]])
  if (!inherits(ssr@fitter ,"function")) stop("code is idiosyncratic for lm fits")
  ps = as.numeric(sapply(ssr, function(x) try(summary(x)$coef[2,4],silent=TRUE)))
}


#plot_mlp = function (ssr, snpMeta, gchr = NULL, geneLocDF=NULL, ps = NULL, pch = 20, cex = 0.5, local = FALSE, plotf=smoothScatter, organism="human") 
#{
#    if (is(ssr@fitter, "GGfitter"))
#        ps = ssr[["pval"]]
#    if (is.null(ps)) 
#        ps = as.numeric(sapply(ssr, function(x) try(summary(x)$coef[2, 
#            4], silent = TRUE)))
#    if (is.null(gchr)) {
#        if (is.null(geneLocDF)) stop("if gchr omitted, please supply geneLocDF")
#	if (!("chr" %in% names(geneLocDF))) stop("geneLocDF does not have chr column.")
#    	gloc = geneLocDF[geneLocDF$gene == ssr@gene, ]
#    	if (nrow(gloc) == 0) 
#        	stop("can't determine gene's chromosome for title from geneLocDF, please supply as arg.")
#    	gchr = gloc[1, "chr"]
#    }
#    bad = which(is.na(ps))
#    if (!local) 
#        XLIM = range(snpMeta[, "pos"])
#    else {
#        snpn = unlist(lapply(ssr, function(x) names(coef(x))[2]))
#        XLIM = range(snpMeta[snpn, "pos"])
#    }
#    if (length(bad)>0) {
#          ssr@locs = ssr@locs[-bad]
#          ps = ps[-bad]
#    }
#    if (!is(snpMeta, "snpMetaWhole")) {
#      plotf(ssr@locs, -log10(ps), xlab = paste("location on chromosome", 
#          chromosome(snpMeta)), ylab = "-log10 p Ho:Bs=0", main = paste(organism, ssr@gene, 
#          "(chr", gchr, ")"), xlim = XLIM, pch = pch, cex = cex)
#      if (!is.null(geneLocDF) & gchr == gsub("chr", "", ssr@chr)) {
#         for (i in 1:nrow(gloc)) {
#             axis(3, at = gloc[i, "beg"], labels = FALSE, col = "green")
#             axis(3, at = gloc[i, "end"], labels = FALSE, col = "red")
#             }
#      }
#    }
#    else if (is(snpMeta, "snpMetaWhole")) {
#      plotf(ssr@locs, -log10(ps), xlab = paste(organism, "chromosome"),
#          ylab = "-log10 p Ho:Bs=0", main = paste(organism, ssr@gene, 
#          "(chr", gchr, ")"), xlim = XLIM, pch = pch, cex = cex, axes=FALSE)
#      axis(2)
#      abline(v=snpMeta@chrbounds, col="gray")
#
#      pts = c(0,snpMeta@chrbounds[-length(snpMeta@chrbounds)]) + diff(c(0,snpMeta@chrbounds))/2
#      #text ( pts, rep(-.1, length(pts) ) , snpMeta@chrlabs )
#      axis(1, at=pts, labels = snpMeta@chrlabs )
#
#      if (!is.null(geneLocDF)) {
#         for (i in 1:nrow(gloc)) {
#             axis(3, at = gloc[i, "beg"], labels = FALSE, col = "green")
#             axis(3, at = gloc[i, "end"], labels = FALSE, col = "red")
#             }
#      }
#    }
#    return(invisible(list(x = ssr@locs, y = -log10(ps))))
#}


plot_mlp = function (ssr, snpMeta, gchr = NULL, geneLocDF = NULL, ps = NULL, 
    pch = 20, cex = 0.6, local = FALSE, plotf = smoothScatter, 
    organism = "human") 
{
    annk = ann = annotation(ssr)
    ann = paste(ann, "db", sep=".")
    annp = paste("package", ann, sep=":")
    require(ann, character.only=TRUE)
    smap = revmap(get(paste(annk, "SYMBOL", sep="")))
    cmap = get(paste(annk, "CHR", sep=""))
    psid = get( ssr@gene, smap )[1]
    if (is.null(gchr)) gchr = chr = get( psid, cmap )[1]
    else chr = gchr
    cname = paste("chr", chr, sep="")
    trans = FALSE
    if (gchr != "all" && cname != snpMeta@chromosome) trans = TRUE
    
    if (is(ssr@fitter, "GGfitter")) 
        ps = ssr[["pval"]]
    if (is.null(ps)) 
        ps = as.numeric(sapply(ssr, function(x) try(summary(x)$coef[2, 
            4], silent = TRUE)))
    if (!is.null(geneLocDF))
        gloc = geneLocDF[geneLocDF$gene == ssr@gene, ]
    bad = which(is.na(ps))
    if (!trans) 
        XLIM = range(snpMeta[, "pos"])
    else {
        XLIM = range(ssr@locs)
    }
    if (length(bad) > 0) {
        ssr@locs = ssr@locs[-bad]
        ps = ps[-bad]
    }
    if (!is(snpMeta, "snpMetaWhole")) {
        plotf(ssr@locs, -log10(ps), xlab = paste("location on chromosome", 
            chromosome(snpMeta)), ylab = "-log10 p Ho:Bs=0", 
            main = paste(organism, ssr@gene, "(chr", gchr, ")"), 
            xlim = XLIM, pch = pch, cex = cex, cex.lab = 1.5, 
            cex.main = 1.5)
        if (!is.null(geneLocDF) & gchr == gsub("chr", "", ssr@chr)) {
            for (i in 1:nrow(gloc)) {
                axis(3, at = gloc[i, "beg"], labels = FALSE, 
                  col = "green")
                axis(3, at = gloc[i, "end"], labels = FALSE, 
                  col = "red")
            }
        }
    }
    else if (is(snpMeta, "snpMetaWhole")) {
        plotf(ssr@locs, -log10(ps), xlab = paste(organism, "chromosome"), 
            ylab = "-log10 p Ho:Bs=0", main = paste(organism, 
                ssr@gene, "(chr", gchr, ")"), xlim = XLIM, pch = pch, 
            cex = cex, axes = FALSE, cex.lab=1.5, cex.main=1.5)
        axis(2)
        abline(v = snpMeta@chrbounds, col = "gray")
        pts = c(0, snpMeta@chrbounds[-length(snpMeta@chrbounds)]) + 
            diff(c(0, snpMeta@chrbounds))/2
        axis(1, at = pts, labels = snpMeta@chrlabs)
        if (!is.null(geneLocDF)) {
            for (i in 1:nrow(gloc)) {
                axis(3, at = gloc[i, "beg"], labels = FALSE, 
                  col = "green")
                axis(3, at = gloc[i, "end"], labels = FALSE, 
                  col = "red")
            }
        }
    }
    return(invisible(list(x = ssr@locs, y = -log10(ps))))
}
