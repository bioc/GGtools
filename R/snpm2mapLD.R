snpm2mapLD = function(x, chrnum, runMAP=TRUE, ...) {
 rsid = colnames(x)
 nsnps = length(rsid)
 sid = rownames(x) # sample
 nsamp = length(sid)
 locs = snpLocs.Hs(chrnum(chrnum), rsid(rsid))["loc",]
 z = as(x, "character")
 if (any(nchar(z) == 0)) {
  warning("for missing genotypes we have substituted A/A")
  z[nchar(z)==0] = "A/A"
 }
 A1 = (apply(z, 1, function(x) sapply(strsplit(x, "/"), "[", 1)))
 A2 = (apply(z, 1, function(x) sapply(strsplit(x, "/"), "[", 2)))
# now individuals are columns, snps are rows
 A1 = as.character(A1)
 A2 = as.character(A2)
 sid = rep(sid, each=nsnps)
 rsid = rep(rsid, nsamp)
 pos = rep(locs, nsamp)
 struc = data.frame(Allele1=A1, Allele2=A2, subjectID=sid, markerID=rsid, position=pos)
 if (runMAP) {
    require(mapLD)
    tmp = mapLD( struc, 4, 3, 1:2, ... )
 }
 else tmp = NA
 list(struc=struc, mapLDans=tmp)
}
