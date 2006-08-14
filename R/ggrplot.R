#ggrplot = function (gge, genesym, snpname = NULL, plat = "hgfocus",
#   psid = NULL) 
#{
#    if (length(snpname) != 1) 
#        stop("must supply single snpname (rsnum)")
#    require(plat, character.only = TRUE)
#    reenv = get(paste(plat, "SYMBOL", sep = ""))
#    rel = as.list(reenv)
#    wh = grep(genesym, unlist(rel))
#    if (length(wh) == 0) stop("could not find requested gene symbol")
#    if (length(wh) > 1 & length(psid) == 0) {
#        warning("more than one probeset for this symbol; using first; set psid to force"); 
#        wh = wh[1]
#    }
#    #gn = get(psid, reenv)
#    if (is.null(psid)) psid = names(rel)[wh]
#    Y = exprs(gge)[psid, ]
#    X = gge@phenoData@pData[[snpname]]
#    plot(X, Y, ylab = genesym, xlab = paste("rare allele count", snpname))
#    mod = lm(Y ~ X)
#    mc = coef(mod)
#    abline(mc[1], mc[2])
#    invisible(list(Y, X, mod))
#}
