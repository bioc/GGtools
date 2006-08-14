#regseq = function(gge, psid, snpinds, snpmeta, ...) {
# getp = function(x){tmp = summary(x)$coef; if(nrow(tmp)==1)NA else tmp[2,4]}
# regs = ggreg(gge, psid, snpinds, ...)
# okreg = regs[ sapply(regs,length) > 1]
# thep = sapply(okreg, getp)
# cpos = snpmeta$pos
# names(cpos) = snpmeta[,1]
# cpploc = thep
# cpploc = cpos[names(thep)]
# list(regs=regs, pvals=thep, locs=cpploc)
#}
