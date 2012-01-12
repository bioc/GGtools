#!/mnt/data/stvjc/ExternalSoft/R-devel-dist/lib64/R/bin/Rscript

suppressPackageStartupMessages(library(GGtools))

cisWriter = function(...) function(ch, stem="ciswrite") {
  dc = function (x) setdiff(as.character(c(1:22, "X", "Y")), x)  # set up a dropChr vector so that only chr x is present
  obj = paste(stem, ch, sep="_")
  fn = paste(obj, ".rda", sep="")
  assign(obj, best.cis.eQTLs(..., dropChr = dc(ch)))
  save(list=obj, file=fn)
  get(obj)
}

ext = function(x) MAFfilter(nsFilter(clipPCs(regressOut(x, ~male), 1:10), var.cutoff=.5), lower=0.05)

library(parallel)
options(mc.cores=13)

demClo = cisWriter("GGdata", ~1, exTransform = ext, unlink.fast=FALSE, geneApply=mclapply)

c2use = commandArgs(TRUE)

set.seed(1234)
if (length(c2use) > 1 | !(c2use %in% as.character(1:22))) stop("unexpected chromosome spec")
demClo(c2use)
