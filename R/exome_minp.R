exome_minpORIG = function( smlSet, fmla, targdir, runname, snpl, mgr=NULL, ... ) {
  if (is.null(mgr) && 
       length(smList(smlSet))>1) stop("requires a chromosome-unified smlSet [length smList == 1]")
  if (!is.null(mgr) && length(mgr@fflist)>1) stop("mgr must have fflist of length 1")
  if (is.null(mgr)) mgr = eqtlTests( smlSet, fmla, targdir=targdir, runname=runname, ... )
  scores = lapply( 1:length(snpl), function(x) mgr[rsid(snpl[[x]]), probeId(names(snpl)[x])] )
  names(scores) = names(snpl)
  ans = sapply( scores, function(x) sapply(x, function(w) apply(w,2,function(z)min(1-pchisq(z, mgr@df)))))
  ans = as.numeric(ans)
  names(ans) = names(snpl)
  ans
}


exome_minp = function (smlSet, fmla, targdir, runname, snpl, feat = NULL,
    mgr = NULL, scoreApply=lapply, ...)
{
    .Deprecated("use ordinary eqtlTests functionality")
    if (is.null(mgr) && length(smList(smlSet)) > 1)
        stop("requires a chromosome-unified smlSet [length smList == 1]")
    if (!is.null(mgr) && length(mgr@fflist) > 1)
        stop("mgr must have fflist of length 1")
    if (is.null(mgr))
        mgr = eqtlTests(smlSet, fmla, targdir = targdir, runname = runname,
            ...)
    if (is.null(feat))
        scores = scoreApply(1:length(snpl), function(x) mgr[rsid(snpl[[x]]),
            probeId(names(snpl)[x])])
    else scores = scoreApply(1:length(snpl), function(x) mgr[rsid(snpl[[x]]),
        probeId(feat)])
    names(scores) = names(snpl)
    ans = sapply(scores, function(x) sapply(x, function(w) apply(w,
        2, function(z) min(1 - pchisq(z, mgr@df)))))
    ans = as.numeric(ans)
    names(ans) = names(snpl)
    ans
}

