cis.FDR.filter.best = function( fn = paste0("comb", 1:22, ".rda"),
    hi.dist = 50000, low.dist = -Inf, hi.maf = .51, low.maf = 0.05,
    fdrOnly = FALSE, applier=lapply ) {
#
# chromsome-specific (or other partition) inclusive runs are filtered
# and post-filter FDR is computed
#
  require(GGtools)
#
# caching
#
  if (!exists(gsub(".rda", "", fn[1]))) objs = lapply(fn, function(x) get(load(x, .GlobalEnv)))
  else objs = lapply(gsub(".rda", "" , fn), get)
  cf = function(x) cisFilter(x, hi.dist = hi.dist, low.dist=low.dist,
        hi.maf=hi.maf, low.maf=low.maf )
  bs = lapply(objs, function(x) {cat("."); bestInStratum(cf(x))})
  bss = lapply(bs, "[[", 1)
  library(parallel)
  ps = applier(1:3, function(x) lapply(objs, function(z) {cat("."); bestInStratum(cf(z), permind=x)}))
  pss = lapply(ps, function(x) lapply(x, function(z) z[[1]]))
  pp = pifdr(unlist(bss), unlist(pss))
  ng = sapply(bs,function(x)length(x[[1]]))
  allg = unlist(lapply(bs, function(x) names(x[[1]])))
  if (fdrOnly) {
     names(pp) = allg
     return(pp)
     } 
  alls = unlist(lapply(bs, "[[", 2))
  allchr = rep(paste0("chr", 1:22), ng)
  ans = data.frame(genes=allg, bestsnp=alls, chr=allchr, fdr=pp)
  ans
}

collectBest = function( fns, targetname="colls", 
   mafs = c(.01, .02, .025, .03333, .05, .075, .1),
   hidists = c(10000, 25000, 50000, 75000, 100000, 250000)) {
nmaf = length(mafs)
ndist = length(hidists)
outli = vector("list", nmaf*ndist) 
opts = as.character(outer(mafs, hidists, paste))
tags = strsplit(opts, " ")
tmp = vector("list", nmaf)
names(tmp) = as.character(mafs)
for (i in 1:nmaf) {
  tmp[[i]] = vector("list", ndist)
  names(tmp[[i]]) = as.character(hidists)
  for (j in 1:ndist) {
     tmp[[i]][[j]] = cis.FDR.filter.best( fns, hi.dist=hidists[j], low.maf=mafs[i] )
     assign(targetname, tmp)
     save(list=targetname, file=paste0(targetname, ".rda"))  # interim
     }
   }
}
