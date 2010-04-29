
sym2id = function (syms, annopkg) 
{
    require(annopkg, character.only = TRUE)
    anbase = gsub(".db", "", annopkg)
    symmap = get(paste(anbase, "SYMBOL", sep = ""))
    unlist(sapply(mget(syms, revmap(symmap), ifnotfound = NA), 
        "[", 1))
}

topFeats = function(probeid=NULL, sym=NULL, rsid=NULL, mgr, ffind, anno, n=10, useSym=TRUE ) {
 if (is.null(probeid) & is.null(sym) & is.null(rsid)) 
       stop("must supply one of probeid, sym, rsid")
 nnn = sum(!is.null(probeid), !is.null(sym), !is.null(rsid))
 if (nnn != 1) stop("only 1 of probeid, sym, rsid must be non-null")
 if (!is.null(sym)) {
  id = sym2id(sym, anno)
  sco = mgr$fflist[[ffind]][, id]
  }
 else if (!is.null(rsid)) sco = mgr$fflist[[ffind]][ rsid, ]
 else if (!is.null(probeid)) sco = mgr$fflist[[ffind]][ ,probeid ]
 inds = order(sco, decreasing=TRUE)[1:n]
 ans = sco[inds]/mgr$shortfac
 if (!is.null(rsid) & useSym) {
  nn = geneSyms(names(ans), anno)
  names(ans) = nn
  }
 ans
}
