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

setClass("hbTestResultsOuter", representation(call="call", hscoreOut="list",
   locs="numeric", formula="ANY", cnum="chrnum", start="numeric", 
   end="numeric", inc="numeric", smlSetName="character",
   blocks="list", loci="list", initSNPs="list"))
setGeneric("pvals", function(x) standardGeneric("pvals"))
setMethod("pvals", "hbTestResultsOuter", function(x) {
 sapply(x@hscoreOut, function(x)x$score.global.p)
})
setMethod("show", "hbTestResultsOuter", function(object) {
 cat("GGtools haplotype block test results.\n")
 cat("Model specified as", object@formula, "on", object@smlSetName, "\n")
 cat("using locations", object@start, "to", object@end, "on chr", object@cnum, "\n")
 cat("Blocks checked in intervals of length", object@inc, "\n")
 cat("There are", length(object@hscoreOut), "test results.\n")
 cat("Minimum global score p:", min(pvals(object)), "\n")
})
 

setGeneric("hbTests", function(fmla, sms, cnum, rsid, rad, ...) 
 standardGeneric("hbTests"))
setMethod("hbTests", c("genesym", "smlSet", "chrnum", "numeric", "numeric"),
  function(fmla, sms, cnum, rsid, rad, ...) {
    theCall = match.call()
    if (length(cnum)!=1) stop("must have length(cnum)==1")
    restr = sms[cnum,]
    ph = as.numeric(exprs(restr[fmla,]))
    sm = smList(restr)[[1]]
    bbHapTests2(ph, cnum, sm, rsid, rad, ...) 
})
#    locs = seq(start, end, inc)
#print(locs)
#    ans = lapply(locs, function(x) {gc(); print(x); 
#       bbHapTests(ph, cnum, sm, x+inc/2, inc/2, doPhase=FALSE, ...)})
#    loci = lapply(ans, function(x) x[[1]]$loci)
#    blocks = lapply(ans, function(x) x[[1]]$blocks)
#    initSNPs = lapply(ans, function(x) x[[1]]$initSNPs)
#    hscoreOut = lapply(ans, function(x) x[[1]]$hapscore)
#    tmp = list(runs=ans, locs=locs)
#    new("hbTestResultsOuter", call=theCall, hscoreOut=hscoreOut,
#     locs=locs, formula=fmla, cnum=cnum, start=start, end=end, inc=inc,
#     smlSetName=deparse(substitute(sms)), blocks=blocks, loci=loci, initSNPs=initSNPs)
#})

  
#bbHapTests2 = function(ph, cnum, sm, rsid, rad=1e5, doPhase=FALSE, ... ) {
#  require(haplo.stats)
#  if (!is(sm, "snp.matrix")) stop("sm must be snp.matrix instance")
#  if (!is(cnum, "chrnum")) stop("cnum must be chrnum instance")
#  slook = snpsNear(rsid, rad, cnum)
#  sm = sm[, rsid(slook)] # has available rsid, must recoup
#  slook = colnames(sm)
#  ldStruc = snpm2mapLD(sm,  cnum, outgraph=tempfile())
#  keep = ldStruc$mapLDans$LDblock$contains
#  inds = sort(as.numeric(unique(unlist(strsplit(keep, ", ")))))
#  fn = slook[inds]
#  smf = sm[, rsid(fn)] 
#  genostruc = snpm2geno(smf)
#  ans = haplo.score(ph, genostruc$geno, locus.label=genostruc$names, ...)
#  schaidRun = list(hapscore=ans, loci=fn, blocks=keep, initSNPs=slook)
##setMethod("invokePhase", c("snp.matrix", "chrnum", "character", "character",
##      "character", "logical"),
##  function(x, cnum, parmstring="",
##     globpname=Sys.getenv("PHASE_LOC"), where2run=".", doParse=TRUE) {
#  if (doPhase) {
#    phobj = invokePhase(smf, cnum, "", Sys.getenv("PHASE_LOC"), ".", TRUE)
#    ph = personalHap(phobj)
#  }
#  else ph = NULL
#  list(schaidRun=schaidRun, phap = ph)
#}



setClass("hbTestResults", representation(hscores="list",
     locs="numeric", chrnum="chrnum", smlSetName="character",
     rsid="ANY", rad="numeric", ldStruc="ANY"))
setGeneric("pvals", function(x)standardGeneric("pvals"))
setMethod("pvals", "hbTestResults", function(x) {
 sapply(x@hscores, function(x)x$score.global.p)
})
setGeneric("locs", function(x) standardGeneric("locs"))
setMethod("locs", "hbTestResults", function(x) x@locs)
setGeneric("hscores", function(x) standardGeneric("hscores"))
setMethod("hscores", "hbTestResults", function(x) x@hscores)
 
setMethod("show", "hbTestResults", function(object){
 cat("GGtools haplotype block test results\n")
 cat("Locations used: ", selectSome(object@locs), "\n")
 cat("Minimum p:", min(pvals(object)), "\n")
})

bbHapTests2 = function (ph, cnum, sm, rsid, rad = 1e+05, ...) 
{
    require(haplo.stats)
    if (!is(sm, "snp.matrix")) 
        stop("sm must be snp.matrix instance")
    if (!is(cnum, "chrnum")) 
        stop("cnum must be chrnum instance")
    slook = snpsNear(rsid, rad, cnum)
    sm = sm[, rsid(slook)]
    slook = colnames(sm)
    ldStruc = snpm2mapLD(sm, cnum, outgraph = tempfile())
    keep = ldStruc$mapLDans$LDblock$contains
    kkk = strsplit(keep, ", ")
    kkl = sapply(kkk, length)
    if (!(any(kkl > 1))) stop("no blocks with length > 1")
    kkk = kkk[kkl>1]
    rs2use = lapply(kkk, function(x) slook[as.numeric(x)])
    ans = list()
    locs = rep(NA, length(rs2use))
    for (j in 1:length(rs2use)) {
      smf = sm[, rsid(rs2use[[j]])]
      genostruc = snpm2geno(smf)
      curloc = mean(snpLocs.Hs(cnum, rsid(rs2use[[j]]))["loc",])
      ans[[j]] = haplo.score(ph, genostruc$geno, 
          locus.label = genostruc$names, ...)
      locs[j] = mnloc=curloc
      }
    new("hbTestResults", hscores=ans, locs=locs, ldStruc=ldStruc,
      smlSetName=deparse(substitute(sm)), rsid=rsid, rad=rad, chrnum=cnum)
}
