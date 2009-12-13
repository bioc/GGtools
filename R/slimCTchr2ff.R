slimCTchr2ff = function (path, chr = 1, dbname = paste("chr", chr, sep = "")) 
{
    require(ff)
    if (!is.numeric(chr)) 
        stop("chr must be numeric index (for Hs, in 1-24)")
    owd = getwd()
    on.exit(setwd(owd))
    tst = try(setwd(path))
    if (inherits(tst, "try-error")) 
        stop("could not change to given folder")
    fn = dir(patt = "\\+.rda$")
    obn = gsub(".rda", "", fn)
    load(fn[1])
    rsn = get(obn[1])$snpnames[[chr]]
    stats = get(obn[1])$chisq[[chr]]
    rm(list = obn[1])
    if (length(obn) > 1) 
        for (i in 2:length(obn)) {
            load(fn[i])
            stats = cbind(stats, get(obn[i])$chisq[[chr]])
            rm(list = obn[i])
        }
    targn = paste(dbname, "chisq",sep="")
    assign(targn, ff(cbind(rsid = rsn, stats), filename=paste(dbname, ".ff", sep=""), dim=c(nrow(stats), ncol(stats)+1),
        dimnames = list(NULL, c("rsid", colnames(stats)))))
    save(list=targn, file=paste(dbname, ".rda", sep=""))
    cat("table written\n")
    NULL
}
