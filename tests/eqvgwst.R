
# april 2011/jan 2012
# compares gwSnpTests to eqtlTests
library(GGtools)
hmceuB36.2021 <- getSS("GGtools", c("20"))
library(illuminaHumanv1.db)
cp = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
hcp = hmceuB36.2021[ probeId(cp), ]
t1 = gwSnpTests(genesym("CPNE1")~male, hcp)
pick = as(t1@.Data[[1]], "data.frame")[22101:22115,]
rsids = rownames(pick)[!is.na(pick[,1])]
csq = pick[rsids,1]
names(csq) = rsids
fi = tempfile()
if (file.exists(fi)) unlink(fi, recursive=TRUE)
t2 = eqtlTests(hcp, ~male, targdir=fi)
sco = t2[rsids,][,1]
unlink(fi, recursive=TRUE)
comp = (sco-trunc(100*csq,0)/100)/sco
(!(max(abs(comp)) > .01))

