 
clipPCs = function (smlSet, inds2drop, center=TRUE)
{
#
# returns smlSet with transformed expressions --
# the principal components in inds2drop are omitted through
# zeroing components of the diagonal component of SVD of t(exprs)
#
    if (!is(smlSet, "smlSet"))
        stop("requires smlSet instance")
    ex = t(exprs(smlSet))
    ex = scale(ex, center=center, scale = FALSE)
    ss = svd(ex)
    d = ss$d
    d[inds2drop] = 0
    recon = t(ss$u %*% diag(d) %*% t(ss$v))
    rownames(recon) = featureNames(smlSet)
    colnames(recon) = sampleNames(smlSet)
    ne = assayDataNew("lockedEnvironment", exprs = recon)
    smlSet@assayData = ne
    smlSet
}

