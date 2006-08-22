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

plot_mlp <- function (ssr, snpMeta, ps=NULL, pch=20, cex=.5) 
{
    if (ssr@fittertok == "fastAGM") ps = ssr[["pval"]]
    else if (ssr@fittertok != "lm") 
        stop("code is idiosyncratic for lm fits")
    if (is.null(ps)) ps = as.numeric(sapply(ssr, function(x) try(summary(x)$coef[2, 
        4], silent = TRUE)))
    if (length(ssr@locs) > 200) plotf = smoothScatter
     else plotf = plot
    data(geneLocs)
    x = geneLocs[geneLocs$gene == ssr@gene, ]
    if (nrow(x) == 0) 
        return(invisible(NULL))
    gchr = x[1,"chr"]
    bad = which(is.na(ps))
    plotf(ssr@locs[-bad], -log10(ps[-bad]), xlab = paste("location on chromosome", 
             chromosome(snpMeta)),
             ylab = "-log10 p Ho:Bs=0", main = paste(ssr@gene,
             "(chr", gchr, ")"), xlim=range(snpMeta[,"pos"]),        
             pch=pch, cex=cex)
    for (i in 1:nrow(x)) {
        axis(3, at = x[i, "beg"], labels = FALSE, col = "green")
        axis(3, at = x[i, "end"], labels = FALSE, col = "red")
    }
    return(invisible(NULL))
}

