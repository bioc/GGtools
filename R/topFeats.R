
sym2id = function (syms, annopkg) 
{
    require(annopkg, character.only = TRUE)
    anbase = gsub(".db", "", annopkg)
    symmap = get(paste(anbase, "SYMBOL", sep = ""))
    unlist(sapply(mget(syms, revmap(symmap), ifnotfound = NA), 
        "[", 1))
}

#setMethod("topFeats", c("rsid"), function(x, ...) {
# .topFeats(probeid=NULL, sym=NULL, rsid=x, ...)
#})
#setMethod("topFeats", c("genesym"), function(x, ...) {
# .topFeats(probeid=NULL, sym=x, rsid=NULL, ...)
#})
#setMethod("topFeats", c("probeId"), function(x, ...) {
# .topFeats(probeid=x, sym=NULL, rsid=NULL, ...)
#})

.topFeats = function(probeid=NULL, sym=NULL, rsid=NULL, mgrOrCTD, ffind=1, anno, n=10, useSym=TRUE, minMAF=0, minGTF=0 ) {
 if (is.null(probeid) & is.null(sym) & is.null(rsid)) 
       stop("must supply one of probeid, sym, rsid")
 if (!is.null(rsid) & (minMAF>0 | minGTF>0)) stop("can't currently use minMAF with rsid selection -- please filter outputs on MAF by hand.")
 if (!is.null(sym)) {
  id = sym2id(sym, anno)
  }
 if (is(mgrOrCTD, "cisTransDirector")) {
  if (minMAF > 0 | minGTF > 0) stop("MAF/GTF filtering not supported for directors yet")
  if (!is.null(sym)) return(sort(mgrOrCTD[, id], decreasing=TRUE)[1:n])
  if (!is.null(rsid)) { tmp = mgrOrCTD[rsid,]; names(tmp)=unlist(sapply(mgrs(mgrOrCTD), function(x)colnames(fffile(x)))); return(sort(tmp, decreasing=TRUE)[1:n] ) }
  if (!is.null(probeid)) return(sort(unlist(mgrOrCTD[,probeid ]), decreasing=TRUE)[1:n])
  }
 nnn = sum(!is.null(probeid), !is.null(sym), !is.null(rsid))
 if (nnn != 1) stop("only 1 of probeid, sym, rsid must be non-null")
 if (!is.null(sym)) {
  sco = fffile(mgrOrCTD)[, id]
  }
 else if (!is.null(rsid)) sco = fffile(mgrOrCTD)[ rsid, ]
 else if (!is.null(probeid)) sco = fffile(mgrOrCTD)[ ,probeid ]
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

.topFeats = function( snpid=NULL, probeid=NULL, mgr, n=10, ffind=1 ) {
  numsupp = sum(c(is.null(snpid), is.null(probeid)))
  if (!(numsupp == 1)) stop("only one of snpid and probeid are allowed")
  curff = mgr@fffile
  on.exit(try(close(mgr@fffile)))
  if (!is.open(curff)) {
     tmp = open(curff)
     if (!is.open(curff)) stop("failed to open ff archive")
     }
  snpn = rownames(curff)
  proben = colnames(curff)
  if (!is.null(snpid)) {
        if (!(snpid %in% snpn)) stop("snpid not found in manager at given ffind")
        tmp = curff[ snpid, ] / mgr@shortfac
        names(tmp) = proben
        }
  else if (!is.null(probeid)) {
        if (!(probeid %in% proben)) stop("probeid not found in manager at given ffind")
        tmp = curff[ , probeid ] / mgr@shortfac
        names(tmp) = snpn
        }
  return(sort(tmp, decreasing=TRUE)[1:n])
}

setMethod("topFeats", c("rsid", "eqtlTestsManager"),
	function(feat, mgr, n=10) {
		.topFeats(snpid=feat, mgr=mgr, n=n, ffind=1)
	})
setMethod("topFeats", c("probeId", "eqtlTestsManager"),
	function(feat, mgr, n=10) {
		.topFeats(probeid=feat, mgr=mgr, n=n, ffind=1)
	})
  
