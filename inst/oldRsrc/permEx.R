
permEx = function(sms) {
#
# permute columns of expression data in an smlSet, ensuring that
# numerical data are assigned to new sample names
#
 ex = exprs(sms)
 nsamp = ncol(ex)
 pinds = sample(1:nsamp, size=nsamp, replace=FALSE)
 ini = colnames(ex)
 pex = ex[,pinds]
 colnames(pex) = ini  # this ensures that sample names are in the
                      # original order
 sms@assayData = assayDataNew("lockedEnvironment", exprs=pex)
 sms
}

