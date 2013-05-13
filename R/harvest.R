cis.FDR.filter.best = function( fn, 
    hi.dist = 50000, low.dist = -Inf, hi.maf = .51, low.maf = 0.05,
    fdrOnly = FALSE, applier=lapply ) {
#
# file-oriented to avoid large internal images
# assumes files are serialized cisRun instances
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
  nperm = length(grep("permScore", names(values(objs[[1]]))))
  cf = function(x) cisFilter(x, hi.dist = hi.dist, low.dist=low.dist,
        hi.maf=hi.maf, low.maf=low.maf )
  chrtags = sapply(objs, function(x) as.character(seqnames(x)[1]))
  bs = lapply(objs, function(x) {cat("."); bestInStratum(cf(x))})
  bss = lapply(bs, "[[", "scores")
  library(parallel)
  ps = applier(1:nperm, function(x) lapply(objs, function(z) {cat("."); bestInStratum(cf(z), permind=x)}))
  pss = lapply(ps, function(x) lapply(x, function(z) z[["scores"]]))
  rawscores = unlist(bss)
  pp = pifdr(rawscores, unlist(pss))
  ng = sapply(bs,function(x)length(x[[1]]))
  allg = unlist(lapply(bs, function(x) names(x[[1]])))
  if (fdrOnly) {
     names(pp) = allg
     return(pp)
     } 
  alls = unlist(lapply(bs, "[[", "scorerids"))  # snpids
  allchr = rep(chrtags, ng)
  ans = data.frame(genes=allg, bestsnp=alls, chr=allchr, fdr=pp, scores=rawscores)
  ans
}

collectBest = function( fns, targetname="harvest", 
   mafs = c(.01, .02, .025, .03333, .05, .075, .1),
   hidists = c(10000, 25000, 50000, 75000, 100000, 250000),
   interimSaves=FALSE) {
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
     if (interimSaves) {
       assign(targetname, tmp)
       save(list=targetname, file=paste0(targetname, ".rda"))  # interim
       }
     }
   }
 tmp
}

  pullHits = function(fns, atts) {
#
# keep at the level of filenames for cisRun instances, and atts
# is the collectBest data frame with optimal selections, filtered to
# desired FDR
#
    tmp = lapply(fns, function(x) get(load(x)))
    kl = lapply(tmp, function(x) paste(names(x), x$snp, sep=":"))
    attk = paste(atts$genes, atts$bestsnp, sep=":")
    tmp = lapply(1:length(tmp), function(x) tmp[[x]][ match( attk, kl[[x]], nomatch=0 ) ])
    curans = do.call(c, lapply(tmp, as, "GRanges"))
#    neword = match( attk, paste(names(curans), curans$snp, sep=":"), nomatch=0)
#   while above seems to work, it appears to be wrong in general
    neword = match( paste(names(curans), curans$snp, sep=":"), attk, nomatch=0) # obtain reordering for atts
    newfdr = atts$fdr[neword]
    curans$fdr = newfdr
    curans$genestart = start(curans)
    curans$geneend = end(curans)
    ranges(curans) = IRanges(curans$snplocs,width=1)
    curans
    }

cis.FDR.filter.SNPcentric = function( fn, 
    hi.dist = 50000, low.dist = -Inf, hi.maf = .51, low.maf = 0.05,
    fdrOnly = FALSE, applier=lapply ) {
#
# file-oriented to avoid large internal images
# assumes files are serialized cisRun instances
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
  nperm = length(grep("permScore", names(values(objs[[1]]))))
  cf = function(x) cisFilter(x, hi.dist = hi.dist, low.dist=low.dist,
        hi.maf=hi.maf, low.maf=low.maf )
  chrtags = sapply(objs, function(x) as.character(seqnames(x)[1]))
  bs = lapply(objs, function(x) {
     x = cf(x)
     tmp = x$score
     names(tmp) = x$probeid
     list(scores=tmp, scorerids=x$snp)
     })
  bss = lapply(bs, "[[", "scores")
  library(parallel)
  pss = applier(1:nperm, function(x) lapply(objs, function(z) {
         values(z)[[paste0("permScore_", x)]]}))
  rawscores = unlist(bss)
  pp = pifdr(rawscores, unlist(pss))
  ng = sapply(bs,function(x)length(x[["scores"]]))
  allg = unlist(lapply(bs, function(x) names(x[["scores"]])))
  if (fdrOnly) {
     names(pp) = allg
     return(pp)
     } 
  alls = unlist(lapply(bs, "[[", "scorerids"))  # snpids
  allchr = rep(chrtags, ng)
  ans = data.frame(genes=allg, snps=alls, chr=allchr, fdr=pp, scores=rawscores)
  ans
}

collectFiltered = function( fns, targetname="harvest", 
   mafs = c(.01, .02, .025, .03333, .05, .075, .1),
   hidists = c(10000, 25000, 50000, 75000, 100000, 250000),
   filterFun = cis.FDR.filter.best,
   interimSaves=FALSE) {
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
     tmp[[i]][[j]] = filterFun( fns, hi.dist=hidists[j], low.maf=mafs[i] )
     if (interimSaves) {
       assign(targetname, tmp)
       save(list=targetname, file=paste0(targetname, ".rda"))  # interim
       }
     }
   }
 tmp
}
