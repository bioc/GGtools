
sym2id = function (syms, annopkg) 
{
    require(annopkg, character.only = TRUE)
    anbase = gsub(".db", "", annopkg)
    symmap = get(paste(anbase, "SYMBOL", sep = ""))
    unlist(sapply(mget(syms, revmap(symmap), ifnotfound = NA), 
        "[", 1))
}

topFeats = function(probeid=NULL, sym=NULL, rsid=NULL, mgrOrCTD, ffind, anno, n=10, useSym=TRUE ) {
 if (is.null(probeid) & is.null(sym) & is.null(rsid)) 
       stop("must supply one of probeid, sym, rsid")
 if (!is.null(sym)) {
  id = sym2id(sym, anno)
  }
 if (is(mgrOrCTD, "cisTransDirector")) {
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
 if (!is.null(rsid) & useSym) {
  nn = geneSyms(names(ans), anno)
  names(ans) = nn
  }
 ans
}