snpm2geno = function (x) 
{
#
# transform a snp.matrix to a "geno" matrix suitable for
# haplo.score, two columns per marker
#
    cnames = colnames(x)
    cc = as(x, "character")
    cc[is.na(cc) | nchar(cc) == 0] = "N/N"
    solWmiss = t(apply(cc, 1, function(x) unlist(strsplit(x, "/"))))
#
# change missing marks in snp.matrix to homozygous common allele
#
    fixN = function(x) {
      if (any(x == "N")) {
         tab = table(x)
         com = names(which.max(tab))
         x[x=="N"] = com[1]
         }
      return(x)
      }
    geno = t(apply(solWmiss, 1, fixN))
    list(geno=geno, names=cnames)
}

bbHapTests = function(ph, cnum, sm, rsid, rad=1e5, doPhase=FALSE, ... ) {
  if (!is(sm, "snp.matrix")) stop("sm must be snp.matrix instance")
  if (!is(cnum, "chrnum")) stop("cnum must be chrnum instance")
  slook = snpsNear(rsid, rad, cnum)
  sm = sm[, rsid(slook)] # has available rsid, must recoup
  slook = colnames(sm)
  ldStruc = snpm2mapLD(sm,  cnum, outgraph=tempfile())
  keep = ldStruc$mapLDans$LDblock$contains
  inds = sort(as.numeric(unique(unlist(strsplit(keep, ", ")))))
  fn = slook[inds]
  smf = sm[, rsid(fn)] 
  genostruc = snpm2geno(smf)
  ans = haplo.score(ph, genostruc$geno, locus.label=genostruc$names, ...)
  schaidRun = list(hapscore=ans, loci=fn, blocks=keep, initSNPs=slook)
#setMethod("invokePhase", c("snp.matrix", "chrnum", "character", "character",
#      "character", "logical"),
#  function(x, cnum, parmstring="",
#     globpname=Sys.getenv("PHASE_LOC"), where2run=".", doParse=TRUE) {
  if (doPhase) {
    phobj = invokePhase(smf, cnum, "", Sys.getenv("PHASE_LOC"), ".", TRUE)
    ph = personalHap(phobj)
  }
  else ph = NULL
  list(schaidRun=schaidRun, phap = ph)
}
  
  
