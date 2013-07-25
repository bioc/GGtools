scoresInStratum = function (m, stratumGetter = names, scoreGetter = function(x) values(x)$score, 
    computeBest = max, scorerName = function(x) values(x)$snp, 
    permind = NULL) 
{
    rng = m
    if (is.null(permind)) 
        sco = scoreGetter(rng)
    else sco = values(rng)[[paste0("permScore_", permind)]]
    sco
}

scoredist = function (fn, hi.dist = 1000, low.dist = -Inf, hi.maf = 0.51, 
    low.maf = 0.01, fdrOnly = FALSE, applier = lapply, perms=FALSE) 
{
    require(GGtools)
    if (!exists(gsub(".rda", "", fn[1]))) 
        objs = lapply(fn, function(x) get(load(x, .GlobalEnv)))
    else objs = lapply(gsub(".rda", "", fn), get)
    nperm = length(grep("permScore", names(values(objs[[1]]))))
    cf = function(x) cisFilter(x, hi.dist = hi.dist, low.dist = low.dist, 
        hi.maf = hi.maf, low.maf = low.maf)
    chrtags = sapply(objs, function(x) as.character(seqnames(x)[1]))
    library(parallel)
    if (!perms) return(scoresInStratum(cf(objs[[1]])))
    ps = applier(1:nperm, function(x) lapply(objs, function(z) {
        cat(".")
        scoresInStratum(cf(z), permind = x)
    }))
    unlist(ps, recursive=TRUE)
}
