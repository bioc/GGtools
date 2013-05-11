

cisFilter = function (rng, low.maf = 0, hi.maf = 0.51, low.dist = 0, hi.dist = 5000 ) {
#
# was initially for mcwAllCis, now can just use cisRun
#
    if (low.dist != 0) stop("low.dist not implemented yet")
    if (!is(rng, "cisRun")) stop("only works for cisRun instance")
    snplocs = rng$snplocs
 # cisRun includes information on radius, and the "range" consists of
 # the inflated radius used for searching
    genestarts = start(rng)+metadata(rng)$radius # shrink back
    geneends = end(rng)-metadata(rng)$radius
    insideCurrentRadius = (snplocs>=(genestarts-hi.dist) &
                           snplocs <= (geneends+hi.dist))
    rng[ which(rng$MAF >= low.maf & rng$MAF < hi.maf & 
        insideCurrentRadius) ]
}

bestInStratum = function (m, stratumGetter = names, scoreGetter = function(x) values(x)$score,
    computeBest = max, scorerName = function(x) values(x)$snp,
    permind = NULL)
{
    rng = m
    if (is.null(permind))
        sco = scoreGetter(rng)
    else sco = values(rng)[[paste0("permScore_", permind)]]
    scorerids = scorerName(rng)
    strat = stratumGetter(rng)
    scoresByStrat = split(sco, strat)
    idsByStrat = split(scorerids, strat)
    bestinds = sapply(scoresByStrat, which.max)
    bestscores = sapply(1:length(bestinds), function(x) scoresByStrat[[x]][bestinds[x]])
    scorers = sapply(1:length(bestinds), function(x) idsByStrat[[x]][bestinds[x]])
    names(bestscores) = names(scoresByStrat)
    list(scores = bestscores, scorerids = as.character(scorers))
}


#bestInStratum = function(m, stratumGetter = names,
#  scoreGetter = function(x) values(x)$score, computeBest = max,
#  scorerName = function(x) values(x)$snp,
#  permind=NULL) {
#  tmp = m
#  if (is.null(permind)) rng = tmp@obs
#  else rng = tmp@perms[[permind]]
#  sco = scoreGetter(rng)
#  scorerids = scorerName(rng)
#  strat = stratumGetter(rng)
#  scoresByStrat = split(sco,strat)
#  idsByStrat = split(scorerids,strat)
#  bestinds = sapply(scoresByStrat, which.max)
#  bestscores = sapply(1:length(bestinds), function(x) scoresByStrat[[x]][bestinds[x]])
#  scorers = sapply(1:length(bestinds), function(x) idsByStrat[[x]][bestinds[x]])
#  names(bestscores)=names(scoresByStrat)
#  list(scores=bestscores, scorerids=as.character(scorers))
#}
  
getFDR = function( files = dir(patt="^pop.*rda$"), filter = cisFilter, 
 reducer = bestInStratum, 
 stratumGetter = names, scoreGetter = function(x) values(x)$score, 
 computeBest = max, ...) {
#
# All.cis oriented approach
#
nfiles = length(files)
perms = obs = vector("list", length(files))
for (i in 1:nfiles) {
  cat(i)
  tmp = get(load(files[i]))
  cf = filter(tmp, ...)
  obs[[i]] = reducer( cf, stratumGetter=stratumGetter,
     scoreGetter=scoreGetter, computeBest = computeBest)
  lobs = lapply(obs, "[[", "scores")
  gnames = lapply(lobs, names)
  nperms = length(tmp@perms)
  perms[[i]] = vector("list", nperms)
  for (j in 1:nperms) {
     obsk = paste(names(cf@obs), cf@obs$snp, sep=":")
     permk = paste(names(cf@perms[[j]]), cf@perms[[j]]$snp, sep=":")
     cf@perms[[j]] = cf@perms[[j]][ match(obsk, permk) ]
     perms[[i]][[j]] = reducer( cf, stratumGetter=stratumGetter,
        scoreGetter=scoreGetter, computeBest = computeBest , permind=j)
     }
  }
  oscores = unlist(lapply(obs, "[[", "scores"))
  pscores = unlist(lapply(perms, function(x) lapply(x, "[[", "scores")), recursive=TRUE)
  fdr = pifdr(oscores, pscores)
# at this point we are not trying to make a comprehensive GRanges
  data.frame(feature=unlist(gnames), fdr=fdr, hitat=unlist(lapply(obs, "[[", "scorerids")))
}

##bigfdr = getFDR()
##save(bigfdr, file="bigfdr.R")
#bigfdr_005_5K = getFDR(hi.dist=5000)
#save(bigfdr_005_5K, file="bigfdr_005_5K.rda")
#bigfdr_005_10K = getFDR(hi.dist=10000)
#save(bigfdr_005_10K, file="bigfdr_005_10K.rda")
#bigfdr_005_50K = getFDR(hi.dist=50000)
#save(bigfdr_005_10K, file="bigfdr_005_10K.rda")
#bigfdr_005_100K = getFDR(hi.dist=100000)
#save(bigfdr_005_100K, file="bigfdr_005_100K.rda")
#bigfdr_025_5K = getFDR( filter = function(x) cisFilter(x, low.maf=.025) )
#save(bigfdr_025_5K, file="bigfdr_025_5K.rda")
#bigfdr_025_10K = getFDR( filter = function(x) cisFilter(x, low.maf=.025, hi.dist=10000) )
#save(bigfdr_025_10K, file="bigfdr_025_10K.rda")
#bigfdr_025_50K = getFDR( filter = function(x) cisFilter(x, low.maf=.025, hi.dist=50000) )
#save(bigfdr_025_50K, file="bigfdr_025_50K.rda")
#bigfdr_025_100K = getFDR( filter = function(x) cisFilter(x, low.maf=.025, hi.dist=100000) )
#save(bigfdr_025_100K, file="bigfdr_025_100K.rda")
#bigfdr_050_5K = getFDR( filter = function(x) cisFilter(x, low.maf=.050) )
#save(bigfdr_050_5K, file="bigfdr_050_5K.rda")
#bigfdr_050_10K = getFDR( filter = function(x) cisFilter(x, low.maf=.050, hi.dist=10000) )
#save(bigfdr_050_10K, file="bigfdr_050_10K.rda")
#bigfdr_050_50K = getFDR( filter = function(x) cisFilter(x, low.maf=.050, hi.dist=50000) )
#save(bigfdr_050_50K, file="bigfdr_050_50K.rda")
#bigfdr_050_100K = getFDR( filter = function(x) cisFilter(x, low.maf=.050, hi.dist=100000) )
#save(bigfdr_050_100K, file="bigfdr_050_100K.rda")
