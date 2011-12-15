
getCSCorrByGene = function(sms, gind, rhs, inicoef=.3) {
 decor = corCompSymm(value=inicoef, form=rhs)
 ex = exprs(sms)[gind,]
 gd = groupedData(as.formula(paste("ex", deparse(rhs))), data=data.frame(ex, pData(sms)))
 decor = Initialize(decor, data=gd)
 tmp = gls(ex~1, cor=decor, data=gd)
 coef(tmp$modelStruct$corStruct, unconstrained=FALSE)
}


