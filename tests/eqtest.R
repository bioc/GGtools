
# this test script exercises eqtlTests and X2chunk extraction of chisq stats
library(GGtools)
try(library(multicore))
applier = lapply
if ("multicore" %in% search()) applier = mclapply
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
library(illuminaHumanv1.db)
gg17 = GGtools:::geneRanges(get("17", revmap(illuminaHumanv1CHR)), "illuminaHumanv1.db")
gnames = gg17[ which(ranges(gg17)$chr17 %in% IRanges(33e6, 35e6)), ]$name
hmlit = hmceuB36.2021[ probeId(gnames), ]
system("rm -rf .abc")
set.seed(1234)
e1 = eqtlTests(hmlit, ~male, targdir=".abc")
data(snpLocs20)
ex2ch33_35live = X2chunk(e1, 1, 33e6, 35e6, snpLocs20, "illuminaHumanv1.db")
load("ex2ch33_35.rda")
all.equal(ex2ch33_35live, ex2ch33_35)
system("rm -rf .abc")
