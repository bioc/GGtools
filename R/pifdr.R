
.pifdrA = function(obs, ps, applier=lapply) {
 nperm=length(ps)/length(obs)
 unlist(applier(obs,
     function(x) (sum(abs(ps)>abs(x))/nperm)/sum(abs(obs)>abs(x))))
}

.pifdrB = function(obs, ps, applier=lapply, use.C=FALSE) {
#
# plug-in FDR -- observed score vector and permuted score vector are inputs
#
        rem = length(ps) %% length(obs)
        if (!(rem==0)) warning("permutation vector length not evenly divisible by observed vector length")
        nperm=length(ps)/length(obs)
#        if (use.C) ( return(.Call("pifdrC", obs, ps, as.integer(nperm))))
        unlist(applier(obs, function(x)
           (sum(abs(ps)>=abs(x))/nperm)/sum(abs(obs)>=abs(x))))
}

pifdr = function(obs, perms, npts=1999, applier=sapply) {
#
# compute plug-in FDR or an approximation to a grid
#   if # observed scores < npts use direct computation
#   otherwise compute for range of points in support of obs
#
  obs = abs(obs)
  perms = abs(perms)
  rem = length(perms) %% length(obs)
  if (!(rem==0)) stop("perms length not integer multiple of obs length")
  nperm = length(perms)/length(obs)
  do.approx = TRUE
  if (length(obs)<=npts) do.approx = FALSE
  if (!do.approx) checkpoints = obs
    else checkpoints = seq(min(obs), max(obs), len=npts)
  calls = applier(checkpoints, function(x) sum(obs >= x))
  false = applier(checkpoints, function(x) sum(perms >= x)/nperm)
  gridFDR = false/calls
  if (!do.approx) return(gridFDR)
  approx(checkpoints, gridFDR, obs, rule=2L)$y
}
  
