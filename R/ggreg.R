ggreg = function (gge, psid, snpinds = NULL, model=c("additive", "general")[1]) 
{
    if (is.null(snpinds)) 
        snpinds = 1:ncol(gge@phenoData@pData)
    out = list()
    j = 1
    for (i in snpinds) {
        if (j%%100 == 0) 
            cat(j)
        sn = gge@phenoData@pData[[i]]
        if (length(unique(sn)) == 1) {
            out[[j]] = NA
            j = j + 1
            next
        }
	if (model == "additive")
        	out[[j]] = lm(exprs(gge)[psid, ] ~ pData(gge)[[i]])
	else if (model == "general")
        	out[[j]] = lm(exprs(gge)[psid, ] ~ factor(pData(gge)[[i]]))
        j = j + 1
    }
    names(out) = names(pData(gge))[snpinds]
    out
}
