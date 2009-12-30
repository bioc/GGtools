library(GGtools)
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
f1 = gwSnpTests(genesym("CPNE1")~male, hmceuB36.2021, chrnum(20))
f2 = gwSnpTests(genesym("PRND")~male, hmceuB36.2021, chrnum(20))
ind = which(f1@.Data[[1]]@test.names == "rs6060535")
sco = f1@.Data[[1]]@chisq[ ind ]
sco2 = f2@.Data[[1]]@chisq[ ind ]
lith = hmceuB36.2021[chrnum(20),]
library(illuminaHumanv1.db)
pid = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
pid2 = get("PRND", revmap(illuminaHumanv1SYMBOL))[1]
pid3 = get("DUSP15", revmap(illuminaHumanv1SYMBOL))[1]
dd = multffCT( list(lith, lith), gs~male, probeId(c(pid,pid2,pid3)))
ddsh = multffCT( list(lith, lith), gs~male, probeId(c(pid,pid2,pid3)), vmode="short", runname="foosh")
getChisq = function(rsid, gene, ctmgr) {
 allrs = lapply(ctmgr$fflist, rownames)
 allg = colnames(ctmgr$fflist[[1]])
 c2use = which(sapply(allrs, function(x) rsid %in% x))
 g2use = which(allg == gene)
 ctmgr$fflist[[c2use]][ rsid, g2use ]
}
# following uses single precision
ee = getChisq("rs6060535", pid, dd)
if (abs(ee/2 - sco) > .0001 ) stop("test failed for scalar gwSnpTests vs multicore multffCT")
ee2 = getChisq("rs6060535", pid2, dd)
if (abs(ee2/2 - sco2) > .0001 ) stop("test failed for scalar gwSnpTests vs multicore multffCT")

# following uses rescaled short int, need bigger tolerance
sumcpne1 = ddsh[rsid("rs6060535"), probeId(pid)][[1]]
if (abs((sumcpne1 - 2*sco)/(2*sco)) > .001) stop("test failed comparing short int based storage to single prec")
