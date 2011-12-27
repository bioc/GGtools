regressOut = function(sms, rhs, ...) {
 mm = model.matrix(rhs, data=pData(sms))
 f = limma::lmFit(exprs(sms), mm, ...)
 r = exprs(sms) - (f$coef %*% t(f$design))
 sms@assayData = assayDataNew("lockedEnvironment", exprs=r)
 sms
}
