#snp.rhs.testsMACH = function(fmla, snp.data=snpdata, data=pData(smlSet), family=      "gaussian") {
#       if (!all.equal(rownames(data), colnames(snp.data))) 
#		stop("rownames of pData must agree with colnames of snp.data")
#       mbase = lm(fmla, data=data)
#       N = apply(snp.data,1,function(x)sum(!is.na(x)))
#       df = rep(1, length(N))
#       cs = sapply(1:nrow(snp.data), function(s) {
#             if (s %% 1000 == 0) cat(s)
#             tmpx <<- snp.data[s,]
#             if (all(tmpx == tmpx[1])) return(NA)
#             fit = try(update(mbase, .~.+tmpx))
#             if (inherits(fit, "try-error")) return(NA)
#             summary(fit)$coef[length(coef(fit)),3]^2  # square the Zstat
#             })
#       new("snp.tests.glm", test.names=rownames(snp.data),
#              chisq=cs, df=as.integer(df), N=as.integer(N))
#}
# 
             
