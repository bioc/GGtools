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
  require(haplo.stats)
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

setClass("hbTestResults", representation(call="call", hscoreOut="list",
   locs="numeric", formula="ANY", cnum="chrnum", start="numeric", 
   end="numeric", inc="numeric", smlSetName="character",
   blocks="list", loci="list", initSNPs="list"))
setGeneric("pvals", function(x) standardGeneric("pvals"))
setMethod("pvals", "hbTestResults", function(x) {
 sapply(x@hscoreOut, function(x)x$score.global.p)
})
setMethod("show", "hbTestResults", function(object) {
 cat("GGtools haplotype block test results.\n")
 cat("Model specified as", object@formula, "on", object@smlSetName, "\n")
 cat("using locations", object@start, "to", object@end, "on chr", object@cnum, "\n")
 cat("Blocks checked in intervals of length", object@inc, "\n")
 cat("There are", length(object@hscoreOut), "test results.\n")
 cat("Minimum global score p:", min(pvals(object)), "\n")
})
 

setGeneric("hbTests", function(fmla, sms, cnum, start, end, inc, ...) 
 standardGeneric("hbTests"))
setMethod("hbTests", c("genesym", "smlSet", "chrnum", "numeric", "numeric", "numeric"),
  function(fmla, sms, cnum, start, end, inc, ...) {
    theCall = match.call()
    if (length(cnum)!=1) stop("must have only a length-1 cnum")
    restr = sms[cnum,]
    ph = as.numeric(exprs(restr[fmla,]))
    sm = smList(restr)[[1]]
    locs = seq(start, end, inc)
print(locs)
    ans = lapply(locs, function(x) {gc(); print(x); 
       bbHapTests(ph, cnum, sm, x+inc/2, inc/2, doPhase=FALSE, ...)})
    loci = lapply(ans, function(x) x[[1]]$loci)
    blocks = lapply(ans, function(x) x[[1]]$blocks)
    initSNPs = lapply(ans, function(x) x[[1]]$initSNPs)
    hscoreOut = lapply(ans, function(x) x[[1]]$hapscore)
    tmp = list(runs=ans, locs=locs)
    new("hbTestResults", call=theCall, hscoreOut=hscoreOut,
     locs=locs, formula=fmla, cnum=cnum, start=start, end=end, inc=inc,
     smlSetName=deparse(substitute(sms)), blocks=blocks, loci=loci, initSNPs=initSNPs)
})

  
