# march 2011 -- test fails as we drop snpLocs metadata.  need to recraft using Herve SNPlocs package
## this test script exercises eqtlTests and X2chunk extraction of chisq stats
#library(GGtools)
#try(library(multicore))
#applier = lapply
#if ("multicore" %in% search()) applier = mclapply
#if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
#library(illuminaHumanv1.db)
#gg17 = GGtools:::geneRanges(get("17", revmap(illuminaHumanv1CHR)), "illuminaHumanv1.db")
#gnames = gg17[ which(ranges(gg17)$chr17 %in% IRanges(33e6, 35e6)), ]$name
#hmlit = hmceuB36.2021[ probeId(gnames), ]
##system("rm -rf .abc")
#set.seed(1234)
#e1 = eqtlTests(hmlit, ~male, targdir=tempdir())
#data(snpLocs20)
#ex2ch33_35live = GGtools:::X2chunk(e1, 1, as.integer(33e6), as.integer(35e6), snpLocs20, "illuminaHumanv1.db")
#load("ex2ch33_35.rda")
#all.equal(ex2ch33_35live, ex2ch33_35)
##system("rm -rf .abc")
#TRUE

# april 2011
library(GGtools)
data(hmceuB36.2021)
library(illuminaHumanv1.db)
cp = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
hcp = hmceuB36.2021[ probeId(cp), ]
hcp = hcp[ chrnum("20"), ]
t1 = gwSnpTests(genesym("CPNE1")~male, hcp, chrnum("20"))
pick = as(t1@.Data[[1]], "data.frame")[22101:22115,]
rsids = rownames(pick)[!is.na(pick[,1])]
csq = pick[rsids,1]
names(csq) = rsids
fi = tempfile()
if (file.exists(fi)) unlink(fi, recursive=TRUE)
t2 = eqtlTests(hcp, ~male, targdir=fi)
sco = t2[rsid(rsids),][[1]][,1]
unlink(fi, recursive=TRUE)
comp = (sco-trunc(100*csq,0)/100)/sco
(!(max(abs(comp)) > .01))


