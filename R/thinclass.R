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
