
setMethod("show", "gwSnpScreenResult", function(object) {
 cat("gwSnpScreenResult with", length(object), "inference data.frames\n")
 cat("gene used: ", object@gene,
    "; expression platform:", object@annotation["exprs"], "\n")
})

#setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setMethod("plot", "gwSnpScreenResult", function(x, y, ...) {
    allp = unlist(lapply(x, "[", , 3))
    snpdata = getSnpData( x@snpLocPackage, x@snpLocExtRef )
    spos = split(as.numeric(snpdata$cumloc), as.numeric(snpdata$chr))
    chrbnd = sapply(spos, max)
#
# following assumes a full human chromosome set ...
#
    allpos = unlist(spos)
    if (length(allp) != length(allpos)) stop(
      "you have used gwSnpScreen without a full set of snp tests; please rerun gwSnpScreen with a chrnum parameter"
      )
    smoothScatter(allpos, -log10(allp), xlab = "genomic position",
        ylab = "-log10 p: chisq1 test", nrp = 200, cex = 0.6,
        pch = 19, main=x@gene)
    mn = apply(cbind(c(0, chrbnd[-length(chrbnd)]), chrbnd), 1, mean)
    clab = c(1:22, "X", "Y")
    abline(v = chrbnd, col = "darkgray")
    axis(3, at = mn, labels = clab, cex = 0.5, las = 2)
})

setMethod("plot", "cwSnpScreenResult", function(x, y, ...) {
    allp = lapply(x, "[", , 3)
    allp = unlist(allp)
    snpAnno = x@annotation["snps"]
#    if (!(exists(snpAnno))) {
#        cat( paste("procuring snpAnno [", snpAnno, "]...\n"))
#        data( list=snpAnno )
#        }
#    gloc = get(snpAnno)
#    loc = gloc$Posi[ which(gloc$Chro == paste("chr", x@chrnum, sep="")) ]
    snpdata = getSnpData( x@snpLocPackage, x@snpLocExtRef )
    loc = snpdata$cumloc[which(as.numeric(snpdata$chr) == x@chrnum)]
    loc = loc - loc[1]
    if (length(x@activeSnpInds) > 0) loc=loc[x@activeSnpInds]
    smoothScatter(loc, -log10(allp), xlab = paste("genomic position on chr",
              x@chrnum),
        ylab = "-log10 p: chisq1 test", nrp = 200, cex = 0.6,
        pch = 19, main=x@gene)
})

setGeneric("topSnps", function(x, ...) standardGeneric("topSnps"))
setMethod("topSnps", "cwSnpScreenResult", function(x, n=10, which="p.1df") {
   x[[1]][ order(x[[1]][,which], decreasing=FALSE), which, drop=FALSE ][1:n,,drop=FALSE]
})

setMethod("topSnps", "gwSnpScreenResult", function(x, n=10, which="p.1df") {
  ts.df = function (w, n = 10, which = "p.1df") 
    {
        w[order(w[, which], decreasing = FALSE), which, 
            drop = FALSE][1:n, , drop = FALSE]
    } 
  lapply(x, ts.df, n=n, which=which)
})



 setMethod("[" , "gwSnpScreenResult", function(x, i, j, ..., drop=FALSE) {
  z = x
  z@.Data = z@.Data[i]
  names(z@.Data) = as.character(i)
  z
 })

setGeneric("getAbsSnpLocs", function(x)standardGeneric("getAbsSnpLocs"))
setMethod("getAbsSnpLocs", c("cwSnpScreenResult"), function(x) {
     met = getSnpData( x@snpLocPackage, x@snpLocExtRef )
     met$cumloc[ met$chr == as(x@chrnum,"numeric") ]
  })
