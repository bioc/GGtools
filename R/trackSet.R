setAs("cwSnpScreenResult", "RangedData", function(from) {
  allp = p.value(from@.Data[[1]]) # , 1) # assume 1df -- must improve
  rs = names(allp)
  locstr = snpLocs.Hs(chrnum(from@chrnum), rsid(rs))
  loc = locstr["loc",]
  locrs = paste("rs", locstr["rsid",], sep="")
  allp = allp[locrs]
  
  require(org.Hs.eg.db, quietly=TRUE)
  rmap = revmap(org.Hs.egSYMBOL)
  ch = paste("chr", from@chrnum, sep="")

  rd = RangedData(IRanges(loc, loc), type = "snpeff", group = "gws",
    score = -log10(allp), space = ch, universe = "hg18")
  allp = rd$score # order probably has changed
  bad = is.na(allp) | !is.finite(allp)
  if (any(bad))
    rd = rd[!bad,]
  rd
})
