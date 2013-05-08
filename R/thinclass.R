setClass("cisRun", contains="GRanges")
setValidity("cisRun", function(object) {
  if (!all(c("MAFlb", "radius", "call", "nperm") %in% names(metadata(object))))
     return("metadata list must include elements MAFlb, radius, and call")
  TRUE
})

convertCis = function(mcw, MAFlb, radius) {
 if (missing(MAFlb)) stop("need MAFlb")
 if (missing(radius)) stop("need radius")
 nperm = length(mcw@perms)
 pnames = paste0("permScore_", 1:nperm)
 ans = mcw@obs
 ansk = paste(names(ans), ans$snp, sep=":")
 for (i in 1:nperm) {
     permk = paste(names(mcw@perms[[i]]), mcw@perms[[i]]$snp, sep=":")
     reord = match(ansk, permk, nomatch=NA)
     if (any(is.na(reord))) stop("permk and ansk don't match")
     values(ans)[[pnames[i]]] = mcw@perms[[i]]$score[reord]
     }
 metadata(ans) = list( MAFlb=MAFlb, radius=radius, call=mcw@theCall, nperm=nperm )
 if (length(metadata(mcw@obs))>0) metadata(ans) = c(metadata(ans), metadata(mcw))
 new("cisRun", ans)
}


hns = function (packnames, chrtag = "chr_22",
   MAFlb, inradius, MAFupdate = function(x,y=1) pmin(x,y))
{
#
# harmonize and sum, for metaanalytic applications with cisRun instances
# the score is simply summed, the MAF is recomputed using MAFupdate, typical
#  approach is to take the minimum MAF among populations
#
    require(GGtools)
    for (i in packnames) require(i, character.only = TRUE)
    objs = lapply(packnames, function(x) get(load(dir(patt = chrtag, 
        system.file("rdas", package = x), full = TRUE))))
    objs = lapply(objs, function(x) {
         if (is(x, "mcwAllCis")) x = convertCis(x, MAFlb=MAFlb, radius=inradius)
         x
         })
    cl = sapply(objs, function(x) is(x, "cisRun"))
    if (!(all(cl))) stop("some inputs not inheriting from cisRun")
    keys = lapply(objs, function(x) paste(names(x), x$snp, 
        sep = ":"))
    okkeys = keys[[1]]
    for (i in 2:length(keys)) okkeys = intersect(okkeys, keys[[i]])
    keepinds = lapply(keys, function(x) match(okkeys,x))
    for (i in 1:length(objs)) objs[[i]] = objs[[i]][keepinds[[i]]]
    newkeys = lapply(objs, function(x) paste(names(x), x$snp, 
        sep = ":"))
    basekeys = newkeys[[1]]
    for (i in 2:length(keys)) {
        if (!(all.equal(newkeys[[i]], basekeys))) stop("keys not matching after filter")
    }
    ans = objs[[1]]
    permcols = grep("^permScore_", names(values(ans)), value=TRUE)
    if (!isTRUE(length(permcols) > 1)) stop("permScore_ entries not found")
    nperms = length(permcols)
    anssco = ans$score
    ansmaf = MAFupdate(ans$MAF)
    # refrain from computing on GRanges but pull the numeric elements
    ansps = vector("list", nperms)
    pnames = vector("character", nperms)
    for (j in 1:nperms) {
        pnames[j] = paste0("permScore_", j)
        ansps[[j]] = values(ans)[[pnames[j]]]
        }
    stopifnot(all(pnames %in% names(values(ans))))
    for (i in 2:length(objs)) {
         anssco = anssco + objs[[i]]$score
         ansmaf = MAFupdate(ansmaf, objs[[i]]$MAF)
         for (j in 1:nperms) {
           ansps[[j]] = ansps[[j]] + values(objs[[i]])[[pnames[j]]]
         }
    }
   ans$score = anssco
   ans$MAF = ansmaf
   for (j in 1:nperms)
        values(ans)[[pnames[j]]] = ansps[[j]]
   ans
}
         
         