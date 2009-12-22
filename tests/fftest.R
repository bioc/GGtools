library(GGtools)
data(hmceuB36.2021)
f1 = gwSnpTests(genesym("CPNE1")~male, hmceuB36.2021, chrnum(20))
ind = which(f1@.Data[[1]]@test.names == "rs6060535")
sco = f1@.Data[[1]]@chisq[ ind ]
lith = hmceuB36.2021[chrnum(20),]
library(illuminaHumanv1.db)
pid = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
dd = multffCT( list(lith, lith), gs~male, probeId(pid))
getChisq = function(rsid, gene, ctmgr) {
 allrs = lapply(ctmgr$fflist, rownames)
 allg = colnames(ctmgr$fflist[[1]])
 c2use = which(sapply(allrs, function(x) rsid %in% x))
 g2use = which(allg == gene)
 ctmgr$fflist[[c2use]][ rsid, g2use ]
}
ee = getChisq("rs6060535", pid, dd)
if (abs(ee/2 - sco) > .0001 ) stop("test failed for scalar gwSnpTests vs multicore multffCT")

