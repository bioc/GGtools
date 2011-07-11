
sym2id = function (syms, annopkg) 
{
    require(annopkg, character.only = TRUE)
    anbase = gsub(".db", "", annopkg)
    symmap = get(paste(anbase, "SYMBOL", sep = ""))
    unlist(sapply(mget(syms, revmap(symmap), ifnotfound = NA), 
        "[", 1))
}

setGeneric("topFeats", function(x, ...)standardGeneric("topFeats"))
setMethod("topFeats", c("rsid"), function(x, ...) {
 .topFeats(probeid=NULL, sym=NULL, rsid=x, ...)
})
setMethod("topFeats", c("genesym"), function(x, ...) {
 .topFeats(probeid=NULL, sym=x, rsid=NULL, ...)
})
setMethod("topFeats", c("probeId"), function(x, ...) {
 .topFeats(probeid=x, sym=NULL, rsid=NULL, ...)
})

.topFeats = function(probeid=NULL, sym=NULL, rsid=NULL, mgrOrCTD, ffind, anno, n=10, useSym=TRUE, minMAF=0, minGTF=0 ) {
 if (is.null(probeid) & is.null(sym) & is.null(rsid)) 
       stop("must supply one of probeid, sym, rsid")
 if (!is.null(rsid) & (minMAF>0 | minGTF>0)) stop("can't currently use minMAF with rsid selection -- please filter outputs on MAF by hand.")
 if (!is.null(sym)) {
  id = sym2id(sym, anno)
  }
 if (is(mgrOrCTD, "cisTransDirector")) {
  if (minMAF > 0 | minGTF > 0) stop("MAF/GTF filtering not supported for directors yet")
  if (!is.null(sym)) return(sort(mgrOrCTD[, id], decreasing=TRUE)[1:n])
  if (!is.null(rsid)) { tmp = mgrOrCTD[rsid,]; names(tmp)=unlist(sapply(mgrs(mgrOrCTD), function(x)colnames(fflist(x)[[1]]))); return(sort(tmp, decreasing=TRUE)[1:n] ) }
  if (!is.null(probeid)) return(sort(unlist(mgrOrCTD[,probeid ]), decreasing=TRUE)[1:n])
  }
 if (is(mgrOrCTD, "multffManager")) {
       fflist = function(x) x$fflist
       shortfac = function(x) x$shortfac
	}
 nnn = sum(!is.null(probeid), !is.null(sym), !is.null(rsid))
 if (nnn != 1) stop("only 1 of probeid, sym, rsid must be non-null")
 if (!is.null(sym)) {
  #id = sym2id(sym, anno)
  #sco = fflist(mgrOrCTD)[[ffind]][, id]
  sco = fflist(mgrOrCTD)[[ffind]][, id]
  }
 else if (!is.null(rsid)) sco = fflist(mgrOrCTD)[[ffind]][ rsid, ]
 else if (!is.null(probeid)) sco = fflist(mgrOrCTD)[[ffind]][ ,probeid ]
 inds = order(sco, decreasing=TRUE)[1:n]
 ans = sco[inds]/shortfac(mgrOrCTD)
# deal with posthoc filtering
 okinds = NULL
     if (minMAF > 0) {
         maf = as.numeric(mgrOrCTD@summaryList[[ffind]][, "MAF"])/mgrOrCTD@shortfac
         okinds = which(maf >= minMAF)
     }
     if (minGTF > 0) {
         mgtf = as.numeric(mgrOrCTD@summaryList[[ffind]][, "mGTF"])/mgrOrCTD@shortfac
         tmp = which(mgtf >= minGTF)
         if (!is.null(okinds)) 
             okinds = intersect(tmp, okinds)
         else okinds = tmp
     }
     if (!is.null(okinds)) {
         okrs = rownames(mgrOrCTD@summaryList[[ffind]])[okinds]
         ans = sort(ans[intersect(okrs, names(ans))],decreasing=TRUE)
     }
#
 if (!is.null(rsid) & useSym) {
  nn = geneSyms(names(ans), anno)
  names(ans) = nn
  }
 ans
}
