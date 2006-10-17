
getfn = function(tok)
 paste("genotypes_chr", tok, "_CEU_r21_nr_fwd.txt.gz", sep="")

alltok = c(1:22,"X", "Y")

getRDAnm = function(tok)
 paste("chr", tok, "GGceuRMA.rda", sep="")


dogz = function(tok, ...) {
 options(verbose=TRUE)
 fn = getfn(tok)
 tmp = HMworkflow(fn, ex, pd, mi, "hgfocus")
 save(tmp, file=getRDAnm(tok), compress=TRUE)
}

library(GGtools)
load("csRMA2.rda")
ex = csRMA2@exprs
data(chr20GGdem)
pd = phenoData(chr20GGdem)
colnames(ex) = rownames(pd@data)
mi = experimentData(chr20GGdem)
for (i in alltok) dogz(i)

