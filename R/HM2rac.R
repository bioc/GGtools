HM2rac = function (fn, comment.char = "", kppref = "NA") 
{
    df = read.table(fn, comment.char = comment.char, h = TRUE)
    allel = as.character(df[, "SNPalleles"])
    kpcol = grep(kppref, names(df))
    df = df[, c(1, kpcol)]
    names(df)[1] = "rsnum"
    BAG = lapply(df, as.character)
    nsnp = nrow(df)
    nsam = ncol(df) - 1
    cmat = matrix(" ", nrow = nsnp, nc = nsam)
    for (i in 1:nsam) cmat[, i] = BAG[[i + 1]]
    names(allel) = rownames(cmat) = BAG[[1]]
    colnames(cmat) = names(df)[-1]
    rr = apply(cmat, 1, countRare)
    rrr = apply(cmat, 1, getRare)
    vdf = data.frame("sampID"=colnames(cmat))
    vmdf = data.frame(labelDescription="sample ID in CEPH system")
    class(vmdf[,1]) = "character"
    vmdf[1,1] = "sample ID in CEPH system"
    rownames(vmdf) = "sampID"
    adf = new("AnnotatedDataFrame", data=vdf, varMetadata=vmdf)
    list(raremat = t(rr), alleles=allel, rareallele = rrr, anno=adf)
}
