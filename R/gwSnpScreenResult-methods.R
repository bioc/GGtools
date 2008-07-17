
setMethod("show", "gwSnpScreenResult", function(object) {
 cat("gwSnpScreenResult with", length(object), "inference data.frames\n")
 cat("gene used: ", object@gene,
    "; expression platform:", object@annotation["exprs"], "\n")
 cat("call was:\n")
 print(object@call)
})

#setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setMethod("plot", "gwSnpScreenResult", function(x, y, doSmooth=TRUE, npts=500, ...) {
#    allp = unlist(lapply(x, "[", , 3))
    allp = unlist(lapply(x, p.value, 1))
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
    if (doSmooth) {
      smoothScatter(allpos, -log10(allp), xlab = "genomic position",
        ylab = "-log10 p: chisq1 test", nrp = npts, cex = 0.6,
        pch = 19, main=x@gene)
        }
    else {
       drop = which(is.na(allp))
       if (length(drop)>0) {
         loc = loc[-drop]
         allp = allp[-drop]
         }
       best = rev(order(-log10(allp)))[1:min(npts,length(loc))]  # inds of top npts 
       plot( allpos[best], -log10(allp)[best], xlab = "genomic position",
         ylab = "-log10 p: chisq1 test", cex = 0.6, pch=19, main=x@gene )
       }
    mn = apply(cbind(c(0, chrbnd[-length(chrbnd)]), chrbnd), 1, mean)
    clab = c(1:22, "X", "Y")
    abline(v = chrbnd, col = "darkgray")
    axis(3, at = mn, labels = clab, cex = 0.5, las = 2)
    tst = try(library(org.Hs.eg.db))
    if (!inherits(tst, "try-error")) {
       rmap = revmap(org.Hs.egSYMBOL)
       egid = get(x@gene, rmap)
       ch = get(egid, org.Hs.egCHR)
       loc = get(egid, org.Hs.egCHRLOC)
       if (ch == "X") ch = 23
        else if (ch == "Y") ch = 24
        else if (ch %in% c("Un", "MT")) ch = 25
        else ch = as.numeric(ch)
       if (ch>1) gpos = chrbnd[ch-1]+abs(loc)
       else gpos = loc
       axis(3, at=gpos, col="red", lwd=2, labels=" ")
       }
})

setMethod("plot", "cwSnpScreenResult", function(x, y, doSmooth=TRUE, npts=500, ...) {
    #allp = lapply(x, "[", , 3)
    allp = p.value(x@.Data[[1]], 1) # assume 1df -- must improve
    #allp = unlist(allp)
    snpAnno = x@annotation["snps"]
#    if (!(exists(snpAnno))) {
#        cat( paste("procuring snpAnno [", snpAnno, "]...\n"))
#        data( list=snpAnno )
#        }
#    gloc = get(snpAnno)
#    loc = gloc$Posi[ which(gloc$Chro == paste("chr", x@chrnum, sep="")) ]
    snpdata = getSnpData( x@snpLocPackage, x@snpLocExtRef )
    offs = getSnpOffsets(x)
    loc = snpdata$cumloc[which(as.numeric(snpdata$chr) == x@chrnum)]
    loc = loc - loc[1] + offs[x@chrnum]
    if (length(x@activeSnpInds) > 0) loc=loc[x@activeSnpInds]
    if (doSmooth) {
      smoothScatter(loc, -log10(allp), xlab = "genomic position",
        ylab = "-log10 p: chisq1 test", nrp = npts, cex = 0.6,
        pch = 19, main=x@gene, xlim=c(0,1.01*max(loc,na.rm=TRUE)))
        }
    else {
         drop = which(is.na(allp))
         if (length(drop)>0) {
            loc = loc[-drop]
            allp = allp[-drop]
            }
       best = rev(order(-log10(allp)))[1:min(npts,length(loc))]  # inds of top npts 
       plot( loc[best], -log10(allp)[best], xlab = "genomic position",
         ylab = "-log10 p: chisq1 test", cex = 0.6, pch=19, main=x@gene ,
         xlim=c(0,1.01*max(loc,na.rm=TRUE)))
       }
       tst = try(library(org.Hs.eg.db))
    if (!inherits(tst, "try-error")) {
       rmap = revmap(org.Hs.egSYMBOL)
       egid = get(x@gene, rmap)
       ch = get(egid, org.Hs.egCHR)
       loc = get(egid, org.Hs.egCHRLOC) # local to chrom!
       axis(3, at=abs(loc), col="red", lwd=2, labels=" ")
       }

#    smoothScatter(loc, -log10(allp), xlab = paste("genomic position on chr",
#              x@chrnum),
#        ylab = "-log10 p: chisq1 test", nrp = 200, cex = 0.6,
#        pch = 19, main=x@gene, xlim=c(0,1.01*max(loc,na.rm=TRUE)))
})

setGeneric("topSnps", function(x, ...) standardGeneric("topSnps"))
setMethod("topSnps", "cwSnpScreenResult", function(x, n=10, which="p.1df") {
   DF = 1
   if (which == "p.2df") DF = 2
   else if (which  != "p.1df") warning("chisq df assumed to be 1, 'which' not recognized")
   pp = p.value(x@.Data[[1]], DF)
   spp = pp[ order(pp, decreasing=FALSE) ]
   df = data.frame(p.val=spp)
   rownames(df) = names(spp)
   df[1:n,,drop=FALSE]
})

setMethod("topSnps", "gwSnpScreenResult", function(x, n=10, which="p.1df") {
  ts.df = function (w, n = 10, which = "p.1df") {
   DF = 1
   if (which == "p.2df") DF = 2
   else if (which  != "p.1df") warning("chisq df assumed to be 1, 'which' not recognized")
   pp = p.value(w, DF)
   spp = pp[ order(pp, decreasing=FALSE) ]
   df = data.frame(p.val=spp)
   rownames(df) = names(spp)
   df[1:n,,drop=FALSE]
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

setMethod("plot", "filteredGwSnpScreenResult", function(x, y, ...) {
 pp = lapply(x@.Data, p.value, 1)
 boxplot(lapply(pp, function(x)-log10(x)), main=x@gene, xlab="chromosome",
   ylab="-log10 p [GLM]")
})

setMethod("plot", "filteredMultiGwSnpScreenResult", function(x, y, ...) {
 stop("please select the desired gene-specific result via [[ and plot directly\n")
})
