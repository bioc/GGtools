
#transTops = function(mgr, chrind=22, n=10) {
#  ranksNvals = ffapply( X=mgr@fflist[[chrind]], AFUN=function(x)c(order(x,decreasing=TRUE)[1:n],
#         sort(x,decreasing=TRUE)[1:n]),
#       MARG=1, N=2*n, RETURN=TRUE, FF_RETURN=FALSE, CFUN="list")
#  if (length(ranksNvals)>1) ranksNvals = do.call("cbind", ranksNvals) else ranksNvals = ranksNvals[[1]]
#  list(ranks=ranksNvals[1:n,,drop=FALSE], vals=ranksNvals[-c(1:n),,drop=FALSE])
#}

