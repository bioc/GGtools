trimSnps = function( listOfSms, rsidlist ) {
  ans = list()
  #cat("making smlSets with common SNPs: \n")
  for (j in 1:length(listOfSms)) {
  #  cat("smlSet", j, "\nchrom ")
    tmp = listOfSms[[j]]
    for (k in 1:length(rsidlist)) {
  #    cat(k)
      tmp@smlEnv$smList[[k]] = listOfSms[[j]]@smlEnv$smList[[k]][, rsidlist[[k]]]
      }
    ans[[j]] = tmp
  #  cat("\n")
    }
  names(ans) = names(listOfSms)
  #cat("done.\n")
  gc()
  ans
}

intersectSnps = function( listOfSms ) {
  nsms = length(listOfSms)
  nchr = length(smList(listOfSms[[1]]))
  chnames = names(smList(listOfSms[[1]]))
  rsidlist = lapply(smList(listOfSms[[1]]), colnames)
  if (length(listOfSms) == 1) return(rsidlist)
  for (j in 2:nsms)
    for (k in 1:length(rsidlist))
      rsidlist[[k]]  = intersect(rsidlist[[k]], colnames(smList(listOfSms[[j]])[[k]]))
  rsidlist
}

makeCommonSNPs = function( listOfSms ) {
  rsidlist = intersectSnps( listOfSms )
  trimSnps( listOfSms, rsidlist )
}


topKfeats = function(mgr, K, fn="inds1.ff", batchsize=200,
   feat = c("score", "ind", "geneind"), ffind=1, ginds ) {
#
# given an eqtlTestsManager instance, presumably generated on one
# chromosome worth of SNP (thus ffind=1 by default), find the
# top K features per SNP in an nsnp x K ff matrix
# 
# to keep things simple, we have separate treatments of values (feat == "score") and
# names of features (feat == "geneind"), where the latter are assumed to be integer
# indices into a vector of probe names
#
# it would of course be very nice to compute the order once and apply it to both
# the scores and the names simultaneously.  i am not sure it will be beneficial
#
     intests = mgr@fffile
     #thevmode = "short"
     #if (feat == "ind" | feat == "geneind") thevmode = "integer"
     if (feat == "score") op = function(x)sort(x, decreasing=TRUE)[1:K]
#     else if (feat == "ind") op = function(x)order(x, decreasing=TRUE)[1:K]
     else if (feat == "geneind") op = function(x)ginds[
                                        order(x,decreasing=TRUE)[1:K]]
     else stop("feat not recognized")
     tmp = ffrowapply(
       t(apply(intests[i1:i2,,drop=FALSE],1,op)),
       X=intests, RETURN=TRUE, RETCOL=K,
       BATCHSIZE=batchsize, VMODE="integer")
     ff(tmp, filename=fn, dim=c(nrow(intests),K), overwrite=TRUE,
       vmode="integer")
       }

.updateKfeats = function( sco1, sco2, ind1, ind2, batchsize=500 ) {
#
# will overwrite sco1, ind1 with the improvements available in sco2, ind2
#
   snchunk = chunk(1, nrow(sco1), by=batchsize)
   K=ncol(sco1)
   for (i in 1:length(snchunk)) {
      i1 = snchunk[[i]][1]
      i2 = snchunk[[i]][2]
      scos = cbind(sco1[i1:i2,,drop=FALSE], sco2[i1:i2,,drop=FALSE])
      ginds = cbind(ind1[i1:i2,,drop=FALSE], ind2[i1:i2,,drop=FALSE])
      rowwiseExtract = function(x,y) t(sapply(1:nrow(x), function(row) x[row,][y[row,]]))
      chind = t(apply(scos, 1, function(x)order(x,decreasing=TRUE)[1:K]))
      sco1[i1:i2,] = rowwiseExtract( scos, chind )
      ind1[i1:i2,] = rowwiseExtract( ginds, chind )
      }
   invisible(NULL)
}


cisZero = function(mgr, snpRanges, geneRanges, radius) {
#
# this function addresses the problem of excluding same-chromosome cis tests
# up to a certain radius around each gene
#
# given a manager, all tests for SNP within radius of each gene are set to zero
#
        dimnt = dimnames(mgr@fffile)
        oks = gsub("rs", "", dimnt[[1]])
        okg = dimnt[[2]]
        srids = values(snpRanges)$RefSNP_id
        if (is.null(srids)) {# perhaps naming convention used
             srids = names(snpRanges)
             if (length(grep("rs", srids))>0) srids = gsub("rs", "", srids)
             }
        if (length(srids) > 0) warning("had to take snp names from names(snpRanges)")
        if (!isTRUE(all(okg %in% names(geneRanges))))  {
            warning("geneRanges does not include names/ranges for all colnames of mgr fffile")
            okg = intersect( okg, names(geneRanges))
            }
        if (!isTRUE(all(oks %in% srids)))  {
            warning("snpRanges does not include names/ranges for all rownames of mgr fffile")
 # match and numeric indexing much more efficient than names for large objects
            oks = match(oks, srids, nomatch=0)
            if (length(oks) == 0) stop("snpRanges does not include any snps in mgr")
            }
 # can't use zero indices in subsetting GRanges...
        ol = findOverlaps(snpRanges[oks[oks>0]], geneRanges[okg] + 
            radius)
        matm = as.matrix(ol)
        if (nrow(matm) > 0) {
            mgr@fffile[matm] = 0
        }
    }

transScores.legacy = function (smpack, snpchr = "chr1", rhs, K = 20, targdirpref = "tsco", 
    geneApply = lapply, chrnames = paste("chr", as.character(1:22), sep=""), 
    geneRanges = NULL, snpRanges = NULL, radius = 2e+06, renameChrs=NULL, 
    probesToKeep=NULL, batchsize=200, genegran=50, shortfac=10, wrapperEndo=NULL,
    geneannopk = "illuminaHumanv1.db", snpannopk = snplocsDefault(),
    gchrpref = "", schrpref="ch", exFilter=function(x)x, 
    smFilter=function(x)x, SSgen=GGBase::getSS)
{
  stopifnot(length(snpchr)==1)

#getCisMap = function( radius=50000, gchr="20",
#  schr="ch20", geneannopk="illuminaHumanv1.db",
#     snpannopk=snplocsDefault(), as.GRangesList=FALSE ) {

#
# objective is a small-footprint accumulation of trans-eQTL tests
#  smpack is a lightweight smlSet package name
#  snpchr is the chromosome for which SNPs will be tested
#  rhs is the right hand side of formula for snp.rhs.tests in snpTests
#  K is the number of best features to be retained (per SNP) as we explore the transcriptome
#  targdirpref is used to define the eqtlTests targdir setting for transient and persistent ff files
#  chrnames is used to enumerate chromosomes for which probes will be identified using the
#    annotation package named in the smlSet -- the "chr" will be stripped
#  geneRanges and snpRanges are used to identify cis tests (relative to radius), which are set to zero
#  renameChrs is available in case the smList elements of the smlSet need to be renamed
#  probesToKeep is used to define a hard restriction on the expression data in the smlSet
#  batchsize is for applications of ffrowapply
#  genegran is for reporting as eqtlTests runs, if options()$verbose is TRUE
#
    #if (length(chrnames) < 2) 
    #    stop("must have length(chrnames) >= 2")
    theCall = match.call()
#
# get an image of the expression+genotype data for SNP on specific chromosome snpchr
#
    sms = SSgen(smpack, snpchr, renameChrs=renameChrs, probesToKeep=probesToKeep, 
       wrapperEndo=wrapperEndo, exFilter=exFilter)
    sms = smFilter(sms)
    if (!is.null(renameChrs)) snpchr=renameChrs
    guniv = featureNames(sms)   # universe of probes
    smanno = gsub(".db", "", annotation(sms))
    require(paste(smanno, ".db", sep = ""), character.only = TRUE)
    clcnames = gsub("chr", "", chrnames)  # typical chrom nomenclature of bioconductor
    pnameList = mget(clcnames, revmap(get(paste(smanno, "CHR", 
        sep = ""))))
 # be sure to use only genes that are on arrays in sms
    pnameList = lapply(pnameList, function(x) intersect(x, guniv))
    names(pnameList) = chrnames  # now the universe of probes is split into chromsomes
    todrop = which(sapply(pnameList, length)==0)
    if (length(todrop)>0) pnameList = pnameList[-todrop]
    genemap = lapply(pnameList, function(x) match(x, guniv))  # numerical indices for probes
    nchr_genes = length(names(pnameList))
    targdir = paste(targdirpref, snpchr, sep="")
#
# start with first element of chrnames vector
#
    inimgr = eqtlTests(sms[probeId(pnameList[[chrnames[1]]]),   # start the sifting through transcriptome
        ], rhs, targdir = targdir, runname = paste("tsc_", chrnames[1],  # testing on genes in chrom 1
        sep = ""), geneApply = geneApply, shortfac=shortfac)
    if (gsub("chr", "", snpchr) == gsub("chr", "", chrnames[1])) {
        mapobj = getCisMap( radius = radius, gchr = paste(gchrpref, chrnames[1], sep=""),
                  schr = paste(schrpref, gsub("chr", "", snpchr), sep=""), geneannopk=geneannopk, snpannopk = snpannopk )
        cisZero(inimgr, mapobj@snplocs, mapobj@generanges, radius=0)   # if SNP are on chrom 1, exclude cis
                             # the gene ranges supplied are already augmented by radius
    }
    topKinds = topKfeats(inimgr, K = K, fn = paste(targdir, "/",  # sort and save
        snpchr, "_tsinds1_1.ff", sep = ""), feat = "geneind", 
        ginds = genemap[[1]], batchsize=batchsize)
    topKscores = topKfeats(inimgr, K = K, fn = paste(targdir, 
        "/", snpchr, "_tssco1_1.ff", sep = ""), feat = "score", 
        ginds = genemap[[1]], batchsize=batchsize)
    kpsnpnames = rownames(inimgr@fffile)
    unlink(filename(inimgr@fffile))
    if (nchr_genes > 1) for (j in 2:nchr_genes) {    # continue sifting through transcriptome
        cat(j)
        gc()
        nxtmgr = eqtlTests(sms[probeId(pnameList[[chrnames[j]]]), 
            ], rhs, targdir = targdir, runname = paste("tsctmp", 
            j, sep = ""), geneApply = geneApply, 
            shortfac=shortfac)
        if (gsub("chr", "", snpchr) == gsub("chr", "", chrnames[j])) {
            mapobj = getCisMap( radius = radius, gchr = paste(gchrpref, chrnames[j], sep=""),
                  schr = paste(schrpref, gsub("chr", "", snpchr), sep=""), geneannopk=geneannopk, snpannopk = snpannopk )
#
# huge bug below, nxtmgr was previously inimgr ... found Aug 21 2012
#
            cisZero(nxtmgr, mapobj@snplocs, mapobj@generanges, radius=0)   # if SNP are on chrom 1, exclude cis
                             # the gene ranges supplied are already augmented by radius
        }
        nxtKinds = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "indscratch.ff", sep = ""), feat = "geneind", ginds = genemap[[j]], 		batchsize=batchsize)
        nxtKscores = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "scoscratch.ff", sep = ""), feat = "score", ginds = genemap[[j]],
                batchsize=batchsize)
        updateKfeats(topKscores, nxtKscores, topKinds, nxtKinds, batchsize=batchsize)  
        unlink(filename(nxtmgr@fffile))   # kill off scratch materials
        unlink(paste(targdir, "indscratch.ff", sep = ""))
        unlink(paste(targdir, "scoscratch.ff", sep = ""))
    }
    baseout = list(scores = topKscores, inds = topKinds, guniv = guniv, K=K, snpchr=snpchr,
        chrnames=chrnames,
	smsanno = annotation(sms),
# previous to 21 aug 2012 snpnames were assigned here from rownames inimgr@fffile; these are now
# memoed above
        snpnames = kpsnpnames, call = theCall, date=date(), shortfac=shortfac)
    new("transManager", base=baseout)
}

updateKfeats = function( sco1, sco2, ind1, ind2, batchsize=200 ) {
#
# will overwrite sco1, ind1 with the improvements available in sco2, ind2
#  this version reduces reliance on ff for sorting tasks
#
   snchunk = chunk(1, nrow(sco1), by=batchsize)
   K=ncol(sco1)
      rsco1 = as.ram(sco1)
      rsco2 = as.ram(sco2)
      rind1 = as.ram(ind1)
      rind2 = as.ram(ind2)
   for (i in 1:length(snchunk)) {
      i1 = snchunk[[i]][1]
      i2 = snchunk[[i]][2]
      scos = cbind(rsco1[i1:i2,,drop=FALSE], rsco2[i1:i2,,drop=FALSE])
      ginds = cbind(rind1[i1:i2,,drop=FALSE], rind2[i1:i2,,drop=FALSE])
      rowwiseExtract = function(x,y) t(sapply(1:nrow(x), function(row) x[row,][y[row,]]))
      chind = t(apply(scos, 1, function(x)order(x,decreasing=TRUE)[1:K]))
      sco1[i1:i2,] = rowwiseExtract( scos, chind )
      ind1[i1:i2,] = rowwiseExtract( ginds, chind )
      }
   rm(rsco1)
   rm(rsco2)
   rm(rind1)
   rm(rind2)
   gc()
   invisible(NULL)
}


bindSnpRanges2mgr = function (mgr, snpRanges, badstart=0L, badend=1e9L) # longer than longest chr
{
#
# this will get a GRanges instance from snpRanges congruent to mgr snp info
# aug 20 2012
# when a SNP is not given location in snpRanges, it will be assumed cis to any gene on same chr
# as it will be given a long interval for findOverlaps
#
    dimnt = dimnames(mgr@fffile)
    snpsInMgr = dimnt[[1]]
    if (!is.null(values(snpRanges)$RefSNP_id)) {
        snpsInMgr = gsub("rs", "", snpsInMgr)
        names(snpRanges) = as.character(values(snpRanges)$RefSNP_id)
    }
    if (is.null(names(snpRanges)))
        stop("need names on snpRanges instances if not a SNPlocs.Hsapiens.* derivative")
    fullsr = IRanges(rep(badstart, length(snpsInMgr)), rep(badend, length(snpsInMgr)))
    sn = seqnames(snpRanges)[1]
    fullsr = GRanges(seqnames = sn, fullsr)
    names(fullsr) = snpsInMgr
    matup = match(names(fullsr), names(snpRanges), nomatch = 0)
    IRanges::end(fullsr[which(matup > 0)]) = IRanges::end(snpRanges[matup[matup>0]])
    IRanges::start(fullsr[which(matup > 0)]) = IRanges::start(snpRanges[matup[matup>0]])
    fullsr
}

bindGeneRanges2mgr = function (mgr, geneRanges, badstart=-6, badend=-5)
{
#
# this will get a GRanges instance from geneRanges congruent to mgr gene info
#
    dimnt = dimnames(mgr@fffile)
    genesInMgr = dimnt[[2]]
    fullgr = IRanges(rep(badstart, length(genesInMgr)), rep(badend, length(genesInMgr)))
    sn = seqnames(geneRanges)[1]
    fullgr = GRanges(seqnames = sn, fullgr)
    names(fullgr) = genesInMgr
    matup = match(names(fullgr), names(geneRanges), nomatch = 0)
    IRanges::end(fullgr[which(matup > 0)]) = IRanges::end(geneRanges[matup[matup>0]])
    IRanges::start(fullgr[which(matup > 0)]) = IRanges::start(geneRanges[matup[matup>0]])
    fullgr
}


cisZero = function (mgr, snpRanges, geneRanges, radius) 
{
#
# use geneRanges+radius to define cis to gene intervals
# find tests of SNPs in these intervals
# set their scores in mgr to zero
#

#getCisMap = function( radius=50000, gchr="20",
#  schr="ch20", geneannopk="illuminaHumanv1.db",
#     snpannopk=snplocsDefault(), as.GRangesList=FALSE ) {

    realSnpRanges = bindSnpRanges2mgr( mgr, snpRanges )
    realGeneRanges = bindGeneRanges2mgr( mgr, geneRanges )
    ol = findOverlaps(realSnpRanges, realGeneRanges +
        radius)
    matm = as.matrix(ol)
    if (nrow(matm) > 0) {
        mgr@fffile[matm] = 0
    }
}

topScores = function(tm) {
 tm@base$scores[,1]/tm@base$shortfac
}
topGenes = function(tm) {
 tm@base$guniv[ tm@base$inds[,1] ]
}
locusNames = function(tm) tm@base$snpnames
geneNames = function(tm) tm@base$guniv
geneIndcol = function(tm, col) tm@base$inds[,col]
nthScores = function(tm, n) {
 tm@base$scores[,n]/tm@base$shortfac
}



setGeneric("transTab", function(x, snps2keep, ...)standardGeneric("transTab"))
setMethod("transTab", c("transManager", "missing"), function(x, snps2keep, ...) {
 .transTab(x@base, NULL, ...)
})
setMethod("transTab", c("transManager", "character"), function(x, snps2keep, ...) {
 .transTab(x@base, snps2keep, ...)
})

.transTab = function( x, snps2keep=NULL, ... ) {
 gannopkname = x$smsanno
 K = x$K
 sids = rep(x$snpnames, each=K )
 gtfs = rep(x$min.gtf, each=K )
 mafs = rep(x$tr.maf, each=K )
 thescos = as.numeric(t(as.ram(x$sco)))/x$shortfac
 theinds = as.numeric(t(as.ram(x$inds)))
 require( gannopkname, character.only=TRUE)
 ganno = gsub(".db", "", gannopkname )
 gloc = sapply(mget(x$guniv, get(paste(ganno, "CHRLOC", sep="")), ifnotfound=NA), "[", 1)
 gchr = sapply(mget(x$guniv, get(paste(ganno, "CHR", sep="")), ifnotfound=NA), "[", 1)
 gsym = sapply(mget(x$guniv, get(paste(ganno, "SYMBOL", sep="")), ifnotfound=NA), "[", 1)
 if (gannopkname != "org.Hs.eg.db") gent = sapply(mget(x$guniv, get(paste(ganno, "ENTREZID", sep="")), ifnotfound=NA), "[", 1)
 else gent = x$guniv
 gn = x$guniv[ theinds ]
 gchr = gchr[ theinds ]
 gsym = gsym[ theinds ]
 gent = gent[ theinds ]
 gloc = gloc[ theinds ]
 okinds = 1:length(sids)
 if (!is.null(snps2keep)) okinds = which(sids %in% snps2keep) 
#
# we need to retrieve the snp locations.  in the trans setting these are only necessary
# when the probe chromosome coincides with snp chromosome, so not always available
# 
 simple = data.frame(snp=sids[okinds], MAF=mafs, GTF=gtfs, chisq=thescos[okinds], probeid=gn[okinds] , probechr=gchr[okinds], snpchr=x$snpchr,
    sym=gsym[okinds], entrez=gent[okinds], geneloc=gloc[okinds], stringsAsFactors=FALSE)
 require(snplocsDefault(), character.only=TRUE)
 stag=x$snpchr[1]  
 stopifnot(all(x$snpchr == stag))
 sl = getSNPlocs(paste0("ch", gsub("chr", "", stag)), as.GRanges=TRUE)
 sln = paste0("rs", sl$RefSNP_id)
 names(sl) = sln
 okn = intersect(simple$snp, sln)
 simple = simple[which(simple$snp %in% okn),]
 starts = start(sl[simple$snp])
 simple$snploc = starts
 data.table(simple)
}


treloc = function (old = "@@REPLACE_BY_RELOCATE@@", new, obj) 
{
# not for export
    tmp = obj
    fref = attr(attributes(tmp@base[["scores"]])[["physical"]], 
        "filename")
    ans = gsub(old, new, fref)
    attr(attributes(tmp@base[["scores"]])[["physical"]], "filename") = ans
    fref = attr(attributes(tmp@base[["inds"]])[["physical"]], 
        "filename")
    ans = gsub(old, new, fref)
    attr(attributes(tmp@base[["inds"]])[["physical"]], "filename") = ans
    obj = tmp
    obj
}


tr1_perm = function () 
{
# for export
    ans = get(load(system.file("transObjs/ptr1.rda", package = "GGtools")))
    treloc("@@REPLACE_BY_RELOCATE@@", system.file(package = "GGtools"), 
        ans)
}

tr1_obs = function () 
{
# for export
    ans = get(load(system.file("transObjs/tr1.rda", package = "GGtools")))
    treloc("@@REPLACE_BY_RELOCATE@@", system.file(package = "GGtools"), 
        ans)
}


mtransScores = function (smpackvec, snpchr = "chr1", rhslist, K = 20, targdirpref = "multtsco", 
    geneApply = lapply, chrnames = paste("chr", as.character(1:22), sep=""), 
    geneRanges = NULL, snpRanges = NULL, radius = 2e+06, renameChrs=NULL,
    batchsize=200, genegran=50, probesToKeep=NULL, shortfac=10, 
    wrapperEndo=NULL, SSgen=GGBase::getSS) 
{
#
# objective is a small-footprint accumulation of trans-eQTL tests 
#      summed over different cohorts
#  smpackvec is a vector of lightweight smlSet package name
#  snpchr is the chromosome for which SNPs will be tested
#  rhslist is the right hand side of formula for snp.rhs.tests in snpTests
#  K is the number of best features to be retained as we explore the transcriptome
#
#
    if (length(chrnames) < 2) 
        stop("must have length(chrnames) >= 2")
    theCall = match.call()
#
# get an image of the expression+genotype data for SNP on specific chromosome snpchr
#
    smsl = lapply(smpackvec, function(x) SSgen(x, snpchr, renameChrs=renameChrs,
       probesToKeep=probesToKeep, wrapperEndo=wrapperEndo))
    if (!is.null(renameChrs)) snpchr=renameChrs
    smsl = makeCommonSNPs(smsl) # could be optional
    names(smsl) = smpackvec
    sms = smsl[[1]]
    guniv = featureNames(sms)   # universe of probes
    smanno = gsub(".db", "", annotation(sms))
    require(paste(smanno, ".db", sep = ""), character.only = TRUE)
    clcnames = gsub("chr", "", chrnames)  # typical chrom nomenclature of bioconductor
    pnameList = mget(clcnames, revmap(get(paste(smanno, "CHR", 
        sep = ""))))
 # be sure to use only genes that are on arrays in sms
    pnameList = lapply(pnameList, function(x) intersect(x, guniv))
    names(pnameList) = chrnames  # now the universe of probes is split into chromsomes
    todrop = which(sapply(pnameList, length)==0)
    if (length(todrop)>0) pnameList = pnameList[-todrop]
    genemap = lapply(pnameList, function(x) match(x, guniv))  # numerical indices for probes
    nchr_genes = length(names(pnameList))
    targdir = paste(targdirpref, snpchr, sep="")
    cursmsl = lapply(smsl, function(x) x[probeId(pnameList[[chrnames[1]]]),])
    inimgr = meqtlTests(cursmsl,   # start the sifting through transcriptome
         rhslist, targdir = targdir, runname = paste("tsc_", chrnames[1],  # testing on genes in chrom 1
        sep = ""), geneApply = geneApply, shortfac=shortfac)
    rm(cursmsl); gc()
    if (gsub("chr", "", snpchr) == gsub("chr", "", chrnames[1])) {
        if (is.null(geneRanges) || is.null(snpRanges)) 
            stop("ranges must be supplied to exclude cis tests")
        cisZero(inimgr, snpRanges, geneRanges, radius)   # if SNP are on chrom 1, exclude cis
    }
    topKinds = topKfeats(inimgr, K = K, fn = paste(targdir, "/",  # sort and save
        snpchr, "_tsinds1_1.ff", sep = ""), feat = "geneind", 
        ginds = genemap[[1]], batchsize=batchsize)
    topKscores = topKfeats(inimgr, K = K, fn = paste(targdir, 
        "/", snpchr, "_tssco1_1.ff", sep = ""), feat = "score", 
        ginds = genemap[[1]], batchsize=batchsize)
    unlink(filename(inimgr@fflist[[1]]))
    for (j in 2:nchr_genes) {    # continue sifting through transcriptome
        cat(j)
        gc()
        cursmsl = lapply(smsl, function(x) x[probeId(pnameList[[chrnames[j]]]),])
        nxtmgr = meqtlTests(cursmsl, rhslist, targdir = targdir, 
            runname = paste("tsctmp", 
            j, sep = ""), geneApply = geneApply, shortfac=shortfac)
        rm(cursmsl); gc()
        if (gsub("chr", "", snpchr) == gsub("chr", "", chrnames[j])) {
            if (is.null(geneRanges) || is.null(snpRanges)) 
                stop("ranges must be supplied to exclude cis tests")
            cisZero(nxtmgr, snpRanges, geneRanges, radius)
            }
        nxtKinds = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "indscratch.ff", sep = ""), feat = "geneind", ginds = genemap[[j]], batchsize=batchsize)
        nxtKscores = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "scoscratch.ff", sep = ""), feat = "score", ginds = genemap[[j]], batchsize=batchsize)
        updateKfeats(topKscores, nxtKscores, topKinds, nxtKinds, batchsize=batchsize)  
        unlink(filename(nxtmgr@fflist[[1]]))   # kill off scratch materials
        unlink(paste(targdir, "indscratch.ff", sep = ""))
        unlink(paste(targdir, "scoscratch.ff", sep = ""))
    }

    baseout = list(scores = topKscores, inds = topKinds, guniv = guniv, K=K, snpchr=snpchr,
        chrnames=chrnames,
	smsanno = annotation(sms),
        snpnames = rownames(inimgr@fffile), call = theCall, date=date(), shortfac=shortfac)
    new("transManager", base=baseout)
}

transScores = function ( tconfig ) {
  snpchr = snpchr(tconfig)
  stopifnot(length(snpchr)==1)
  smpack = smpack(tconfig)
  rhs = rhs(tconfig)
  K = gbufsize(tconfig)
  targdirpref = folderStem(tconfig)
  geneApply = geneApply(tconfig)
  chrnames = chrnames(tconfig)
  renameChrs = paste0(smchrpref(tconfig), chrnames)
  radius = radius(tconfig)
  geneRanges = snpRanges = NULL
  batchsize = batchsize(tconfig)
  shortfac = shortfac(tconfig)
  probesToKeep = wrapperEndo = NULL
  geneannopk = geneannopk(tconfig)
  snpannopk = snpannopk(tconfig)
  gchrpref = gchrpref(tconfig)
  schrpref = schrpref(tconfig)
  exFilter = exFilter(tconfig)
  smFilter = smFilter(tconfig)
  SSgen = SSgen(tconfig)
  
    theCall = match.call()
#
# get an image of the expression+genotype data for SNP on specific chromosome snpchr
#
    sms = SSgen(smpack, snpchr, renameChrs=snpchr, probesToKeep=probesToKeep, 
       wrapperEndo=wrapperEndo, exFilter=exFilter)
    sms = smFilter(sms)
#    if (!is.null(renameChrs)) snpchr=renameChrs
    guniv = featureNames(sms)   # universe of probes
    smanno = gsub(".db", "", annotation(sms))
    require(paste(smanno, ".db", sep = ""), character.only = TRUE)
    clcnames = gsub("chr", "", chrnames)  # typical chrom nomenclature of bioconductor
# NOOP if no prefix to begin with
    pnameList = mget(clcnames, revmap(get(paste(smanno, "CHR", 
        sep = ""))))
 # be sure to use only genes that are on arrays in sms
    pnameList = lapply(pnameList, function(x) intersect(x, guniv))
    names(pnameList) = chrnames  # now the universe of probes is split into chromsomes
    todrop = which(sapply(pnameList, length)==0)
    if (length(todrop)>0) pnameList = pnameList[-todrop]
    genemap = lapply(pnameList, function(x) match(x, guniv))  # numerical indices for probes
    nchr_genes = length(names(pnameList))
    targdir = paste(targdirpref, snpchr, sep="")
#
# start with first element of chrnames vector
#
    inimgr = eqtlTests(sms[probeId(pnameList[[chrnames[1]]]),   # start the sifting through transcriptome
        ], rhs, targdir = targdir, runname = paste("tsc_", chrnames[1],  # testing on genes in chrom 1
        sep = ""), geneApply = geneApply, shortfac=shortfac)
    if (gsub("chr", "", snpchr) == gsub("chr", "", chrnames[1])) {
        mapobj = getCisMap( radius = radius, gchr = paste(gchrpref, chrnames[1], sep=""),
                  schr = paste(schrpref, gsub("chr", "", snpchr), sep=""), geneannopk=geneannopk, snpannopk = snpannopk )
        cisZero(inimgr, mapobj@snplocs, mapobj@generanges, radius=0)   # if SNP are on chrom 1, exclude cis
                             # the gene ranges supplied are already augmented by radius
    }
    topKinds = topKfeats(inimgr, K = K, fn = paste(targdir, "/",  # sort and save
        snpchr, "_tsinds1_1.ff", sep = ""), feat = "geneind", 
        ginds = genemap[[1]], batchsize=batchsize)
    topKscores = topKfeats(inimgr, K = K, fn = paste(targdir, 
        "/", snpchr, "_tssco1_1.ff", sep = ""), feat = "score", 
        ginds = genemap[[1]], batchsize=batchsize)
    kpsnpnames = rownames(inimgr@fffile)
    unlink(filename(inimgr@fffile))
    if (nchr_genes > 1) for (j in 2:nchr_genes) {    # continue sifting through transcriptome
        cat(j)
        gc()
        nxtmgr = eqtlTests(sms[probeId(pnameList[[chrnames[j]]]), 
            ], rhs, targdir = targdir, runname = paste("tsctmp", 
            j, sep = ""), geneApply = geneApply, 
            shortfac=shortfac)
        if (gsub("chr", "", snpchr) == gsub("chr", "", chrnames[j])) {
            mapobj = getCisMap( radius = radius, gchr = paste(gchrpref, chrnames[j], sep=""),
                  schr = paste(schrpref, gsub("chr", "", snpchr), sep=""), geneannopk=geneannopk, snpannopk = snpannopk )
#
# huge bug below, nxtmgr was previously inimgr ... found Aug 21 2012
#
            cisZero(nxtmgr, mapobj@snplocs, mapobj@generanges, radius=0)   # if SNP are on chrom 1, exclude cis
                             # the gene ranges supplied are already augmented by radius
        }
        nxtKinds = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "indscratch.ff", sep = ""), feat = "geneind", ginds = genemap[[j]], 		batchsize=batchsize)
        nxtKscores = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "scoscratch.ff", sep = ""), feat = "score", ginds = genemap[[j]],
                batchsize=batchsize)
        updateKfeats(topKscores, nxtKscores, topKinds, nxtKinds, batchsize=batchsize)  
        unlink(filename(nxtmgr@fffile))   # kill off scratch materials
        unlink(paste(targdir, "indscratch.ff", sep = ""))
        unlink(paste(targdir, "scoscratch.ff", sep = ""))
    }
    ss.gtf = GGBase:::smlSummary(sms)[[1]]  # only one chrom
    ss.maf = ss.gtf[,"MAF"]
    names(ss.maf) = rownames(ss.gtf)
    mingtfs = apply(ss.gtf, 1, function(z) min(z[c("P.AA", "P.AB", "P.BB")]))
    names(mingtfs) = rownames(ss.gtf)
    if (length(setdiff(kpsnpnames, names(mingtfs)))>0) warning("some kept snp  had no GTF")
    min.gtf = mingtfs[kpsnpnames]
    tr.maf = ss.maf[kpsnpnames]

    baseout = list(scores = topKscores, inds = topKinds, guniv = guniv, K=K, snpchr=snpchr,
        chrnames=chrnames,
	smsanno = annotation(sms),
# previous to 21 aug 2012 snpnames were assigned here from rownames inimgr@fffile; these are now
# memoed above
        snpnames = kpsnpnames, call = theCall, date=date(), shortfac=shortfac, min.gtf=min.gtf, tr.maf=tr.maf,
        config = tconfig)
    new("transManager", base=baseout)  # note includes references to ff in baseout$scores
}


cleanup_transff = function(x) {
 fn = attr(attr(x@base$scores, "physical"), "filename")
 comps = strsplit(fn, "/")[[1]]
 nel = length(comps)
 unlink(comps[nel-1], recursive=TRUE)
}


transeqByCluster = function( cl, snpchrs=c("chr21", "chr22"), exchrs=1:22, baseconf, targname="transrun_", nperm=1, inseed=1234, ... ) {
 RNGkind("L'Ecuyer-CMRG")
 set.seed(inseed)
 baseconf <<- baseconf
 targname <<- targname
 clusterExport(cl, "baseconf")
 clusterExport(cl, "targname")
 assign(cfn <- paste0(targname, "baseconf"), baseconf)
 save(list=cfn, file=paste0(cfn, ".rda"))
 ttabs = clusterApplyLB( cl, snpchrs, function(x) {
    gc() 
    snpchr(baseconf) = x # only manipulation apart from init prior to cal
    tab <- transTab( tmp <- transScores( baseconf ) ) 
    cleanup_transff(tmp) 
    pscolist = vector("list", nperm)
    for (k in 1:nperm) {
       smFilter(baseconf) = permEx
       pscolist[[k]] = transTab( tmp <- transScores( baseconf ) )$chisq
       cleanup_transff(tmp) 
       }
    pnames = paste0("permScore_", 1:nperm)
    for (k in 1:nperm) tab[[pnames[k]]] = pscolist[[k]]
    obn = paste0(targname, x)
    assign(obn, tab)
    save(list=obn, file=paste0(obn, ".rda"))  # one data table per chrom
 } )
 invisible(NULL)
}

