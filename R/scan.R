eqtlscan = function( 
    infmla=`GI_23397697-A`~1, smpack="GGdata", chrtag="20", ...) {
#
# decode formula 
#
    respObj = infmla[[2]] # we know infmla is a formula, infmla[[2]] is dep var
#
# at this point we have the featureName that we need
#
    pname = as.character(respObj)
    infmla[[2]] = as.name(pname)  # replace the dependent variable spec in fmla
    ex = NULL  # scotch global variable NOTE
    data(eset, package=smpack) # defines ex
    gtdat = get(load(system.file(paste("parts/", chrtag, ".rda", sep="", collapse=""),
      package=smpack)))
    assign(pname, exprs(ex)[pname,]) # expression vector
    alld = data.frame(get(pname), pData(ex))
    names(alld)[1] = pname
    snp.rhs.tests(infmla, family="gaussian",
        snp.data=gtdat, data=alld, uncertain=TRUE)
}
