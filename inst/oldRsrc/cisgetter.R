
getSNPgr = function(genome, chr) {
 if (!(chr %in% as.character(1:22))) stop("chr must be in as.character(1:22)")
 if (genome == "hg18") spack = "SNPlocs.Hsapiens.dbSNP.20090506"
 else if (genome == "hg19") spack = "SNPlocs.Hsapiens.dbSNP.20100427"
 else stop("only using hg18 and hg19 to select SNPlocs")
 require(spack, character.only=TRUE)
 if (packageVersion(spack) < package_version("0.99.6")) stop(
        "please install SNPlocs package with version at least 0.99.6 [these have GRanges conversion built in]")
 seqn = do.call(":::", list(spack, "SEQNAMES"))  #
# here we will assume that prefix of seqn is either "ch" or "chr"
 nseqn = gsub("ch", "", gsub("chr", "", seqn))
 seqnToUseInd = which(nseqn == chr)
 if (length(seqnToUseInd)==0) stop(paste("no match of chr submitted in call to SEQNAMES of", spack))
 seqnToUse = seqn[seqnToUseInd]
 getter = do.call("::", list(spack, "getSNPlocs"))  # for multiple SNPlocs on searchlist
 snpgr = getter( seqnToUse, as.GRanges=TRUE )  # LAST USE of seqnToUse, reset to chr prefix below
}

best.cis.eQTLs = function(smpack, fmla, cisRadius=50000, genome="hg19", 
   exTransform = function(x)x, folderstem="cisScratch", nperm=2, 
   geneApply=lapply, chromApply=lapply, cleanChrn = function(x) gsub("chr", "", x),
   additionalSNPGR=NULL, useTxDb=FALSE, verbose=TRUE,
   dropChr = c("X", "Y", "M"), unlink.fast=TRUE) {
#
# this program takes an smlSet package and identifies for each probe
# the best "cis" associated SNP.  probes and SNPs will be filtered using exTransform
#
# the outputs are to be incrementable
#
  if (!identical(chromApply, lapply) && unlink.fast==TRUE) {
      warning("you seem to have non-sequential iteration over chromosomes, setting unlink.fast to FALSE")
      unlink.fast = FALSE
      }
#
# phase 1 obtain appropriate chromosome labels for the genome-wide run
#
  smparts = dir(system.file("parts", package=smpack))
  chrn = gsub(".rda", "", smparts)
  cchrn = cleanChrn(chrn)
  bad = which(cchrn %in% dropChr)
  if (length(bad) != 0) {
     cchrn = cchrn[-bad]
     chrn = chrn[-bad]
     }
  nchr = length(chrn)
  if (verbose) cat("will run on chr\n")
  print(chrn)
#
# iterate over chromosomes
#
  ans = chromApply(1:nchr, function(i) {
  #for (i in 1:nchr) {
    if (verbose) cat("chr", chrn[i], "... getSS...")
    cursms = getSS(smpack, chrn[i])
    if (!require(annotation(cursms), character.only=TRUE)) 
      stop(paste("got annotation(sms) == ",
          annotation(cursms), "but need a .db annotation package\n",
	  " for the expression component of smpack"))
    if (verbose) cat("exTransform ...")
    # this exTrans operation could be done once at top, but should be fast...
##
## filter and restrict probes
##
    cursms = exTransform(cursms) # should work on full expression component,
			# for best heterogenity reduction, var filtering...
    cursms = restrictProbesToChrom( cursms, cchrn[i] ) # should be late
    if (verbose) cat("gene2snp map ...\n")
##
## computations for cis proximity
##
    ge = geneRanges( annotation(cursms), paste("chr", cchrn[i], sep=""),
       is.annopkg=TRUE )
    sr = getSNPgr( genome, cchrn[i] )
    g2sl = getGene2SnpList( cursms, cchrn[i], genome=genome, radius=cisRadius,
	additionalSNPGR=additionalSNPGR, useTxDb=useTxDb)
##
## direct testing and permutations saved in FDRtab output list
##
    tmp = genewiseFDRtab(cursms, fmla,
        geneApply=geneApply, chromApply=lapply, # at this point you are just running one chromosome
        folderstem=folderstem, nperm=nperm, 
        geneExtents=ge, snpRanges=sr, force.locations=TRUE,
        gene2snpList=g2sl)
##
## now clean up
##
    dirs = paste("p", 1:nperm, folderstem, sep="_")
    dirs = c(folderstem, dirs)
    if (identical(chromApply, lapply) && unlink.fast==TRUE) unlink(dirs, recursive=TRUE, force=TRUE)
    gc()
    tmp
    })
   tmp = do.call(c, ans)  # this binds lists and converts structure to eqtlFDRSummary
   tmp@theCall = match.call()
   tmp@sess = sessionInfo()
   tmp@genome=genome
   tmp@cisRadius=cisRadius
   tmp@nperm=nperm
   tmp@exTransform = exTransform
   if (!unlink.fast) {
    dirs = paste("p", 1:nperm, folderstem, sep="_")
    dirs = c(folderstem, dirs)
    unlink(dirs, recursive=TRUE, force=TRUE)
    }
   tmp
}

all.cis.eQTLs = function(smpack, fmla, cisRadius=50000, genome="hg19", 
   exTransform = function(x)x,
   folderstem="cisScratch", 
   geneApply=lapply, chromApply=lapply,
   cleanChrn = function(x) gsub("chr", "", x),
   additionalSNPGR=NULL, useTxDb=FALSE, verbose=TRUE,
   dropChr = c("X", "Y", "M"), unlink.fast=TRUE) {
  if (!identical(chromApply, lapply) && unlink.fast==TRUE) {
      warning("you seem to have non-sequential iteration over chromosomes, setting unlink.fast to FALSE")
      unlink.fast = FALSE
      }
  smparts = dir(system.file("parts", package=smpack))
  chrn = gsub(".rda", "", smparts)
  cchrn = cleanChrn(chrn)
  bad = which(cchrn %in% dropChr)
  if (length(bad) != 0) {
     cchrn = cchrn[-bad]
     chrn = chrn[-bad]
     }
  nchr = length(chrn)
  if (verbose) cat("will run on chr\n")
  print(chrn)
  g2smaps = list()
  ans = chromApply(1:nchr, function(i) {
  #for (i in 1:nchr) {
    if (verbose) cat("chr", chrn[i], "... getSS...")
    cursms = getSS(smpack, chrn[i])
    if (!require(annotation(cursms), character.only=TRUE)) 
      stop(paste("got annotation(sms) == ",
          annotation(cursms), "but need a .db annotation package\n",
	  " for the expression component of smpack"))
    if (verbose) cat("exTransform ...")
    # this exTrans operation could be done once at top, but should be fast...
    cursms = exTransform(cursms) # should work on full expression component,
			# for best heterogenity reduction, var filtering...
    cursms = restrictProbesToChrom( cursms, cchrn[i] ) # should be late
    if (verbose) cat("gene2snp map ...\n")
    ge = geneRanges( annotation(cursms), paste("chr", cchrn[i], sep=""),
       is.annopkg=TRUE )
    sr = getSNPgr( genome, cchrn[i] )
    g2sl = getGene2SnpList( cursms, cchrn[i], genome=genome, radius=cisRadius,
	additionalSNPGR=additionalSNPGR, useTxDb=useTxDb)
    tmp = eqtlTests(cursms, fmla,
        geneApply=geneApply, chromApply=lapply, # at this point you are just running one chromosome
        targdir=folderstem)
    gn = probesManaged(tmp, 1)
    okrs = snpsManaged(tmp, 1)
    g2sl = g2sl[intersect(names(g2sl), gn)]
    g2sl = lapply( g2sl, function(x) intersect(x, okrs))
    tmp = lapply(names(g2sl), function(x) tmp[rsid(g2sl[[x]]), probeId(x) ])
    unlink(folderstem, recursive=TRUE)
    tmp
    })
   ans
}
