
setMethod("trackSet", "cwSnpScreenResult", function(object, ...) {
    allp = p.value(object@.Data[[1]]) # , 1) # assume 1df -- must improve
    rs = names(allp)
    locstr = snpLocs.Hs(chrnum(object@chrnum), rsid(rs))
    loc = locstr["loc",]
    locrs = paste("rs", locstr["rsid",], sep="")
    allp = allp[locrs]
 
    require(org.Hs.eg.db, quietly=TRUE)
    rmap = revmap(org.Hs.egSYMBOL)
    ch = paste("chr", object@chrnum, sep="")

    fdata = data.frame(chrom=ch, start=loc, end=loc, strand="+", phase=NA,
       type="snpeff", group="gws", score=-log10(allp))
    bad = which(is.na(allp) | !is.finite(allp))
    if (length(bad) > 0) {
     adata = assayDataNew("lockedEnvironment", dataVals=matrix(-log10(allp[-bad]),nc=1))
     fd = new("AnnotatedDataFrame", data=fdata[-bad,])
    } else {
     adata = assayDataNew("lockedEnvironment", dataVals=matrix(-log10(allp),nc=1))
     fd = new("AnnotatedDataFrame", data=fdata)
    }
    require(rtracklayer, quietly=TRUE)
    new("trackSet", genome="Hs", assayData=adata, featureData=fd)
})

