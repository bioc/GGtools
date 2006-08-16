
fastAGM = function(snpm, exprv) {
 if (any(is.na(snpm))) stop("some missing genotypes; should not get to this point")
 ans = .C("mreg_engine", as.integer(nrow(snpm)),
    as.integer(ncol(snpm)), as.double(exprv), as.double(snpm),
    b=double(nrow(snpm)), se=double(nrow(snpm)))
 trat = ans$b/ans$se
 pval = pmin(1,2*pt(-abs(trat),df=length(exprv)-2))
 names(ans$b) = names(ans$se) = names(trat)= names(pval) = rownames(snpm)
 list(b=ans$b, se=ans$se, trat=trat, pval=pval)
}

