
setGeneric("filterSnpTests", function(x, n) standardGeneric(
 "filterSnpTests"))

topSingSnps = function(x, n=250, df=1) {
  pp = p.value(x, df)
  opp = order(pp)
  x[opp[1:n]]
}

filterGWS = function(x, ...) {
  x@.Data = lapply(x@.Data, topSingSnps, ...)
  new("filteredGwSnpScreenResult", x)
}


setMethod("filterSnpTests",
  "multiGwSnpScreenResult", function(x, n) {
    tmp = lapply( x@.Data, filterGWS, n )
    names(tmp) = geneIds(x@geneset)
    new("filteredMultiGwSnpScreenResult", geneset=x@geneset,
      call=x@call, tmp)
})

setMethod("filterSnpTests",
  "gwSnpScreenResult", function(x, n) {
   filterGWS(x, n)
})


setMethod("plot", "filteredGwSnpScreenResult", function(x, y, ...) {
 pp = lapply(x@.Data, p.value, 1)
 boxplot(lapply(pp, function(x)-log10(x)), main=x@gene, xlab="chromosome",
   ylab="-log10 p [GLM]")
 xx = try(require(org.Hs.eg.db, quietly=TRUE))
 if (!inherits(xx, "try-error")) {
  if (is(x@gene, "genesym")) rmap = revmap(org.Hs.egSYMBOL)
 else if (is(x@gene, "probeId"))  {
           require(x@annotation, character.only=TRUE, quietly=TRUE)
           rmap = get(paste(gsub(".db", "", x@annotation, "ENTREZID", sep="")))
           }

  else {
    warning("x@gene is neither symbol nor probeID, we do not plot the location.")
    return(invisible(NULL))
    }
    egid = get(x@gene, rmap)
    ch = try(get(egid, org.Hs.egCHR))
    if (!inherits(ch, "try-error")) {
       if (ch == "X") ch = 23
       else if (ch == "Y") ch = 24
       axis(3, at=as.numeric(ch), col="red", labels=" ")
       }
    }
})

setMethod("plot", "filteredMultiGwSnpScreenResult", function(x, y, ...) {
 stop("please select the desired gene-specific result via [[ and plot directly\n")
})
