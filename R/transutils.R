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
     intests = mgr@fflist[[ffind]]
     #thevmode = "short"
     #if (feat == "ind" | feat == "geneind") thevmode = "integer"
     if (feat == "score") op = function(x)sort(x, decreasing=TRUE)[1:K]
#     else if (feat == "ind") op = function(x)order(x, decreasing=TRUE)[1:K]
     else if (feat == "geneind") op = function(x)ginds[
                                        order(x,decreasing=TRUE)[1:K]]
     else stop("feat not recognized")
     tmp = ffrowapply(
       t(apply(intests[i1:i2,],1,op)),
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
      scos = cbind(sco1[i1:i2,], sco2[i1:i2,])
      ginds = cbind(ind1[i1:i2,], ind2[i1:i2,])
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
        dimnt = dimnames(mgr@fflist[[1]])
        oks = gsub("rs", "", dimnt[[1]])
        okg = dimnt[[2]]
        srids = values(snpRanges)$RefSNP_id
        if (is.null(srids)) {# perhaps naming convention used
             srids = names(snpRanges)
             if (length(grep("rs", srids))>0) srids = gsub("rs", "", srids)
             }
        if (length(srids) > 0) warning("had to take snp names from names(snpRanges)")
        if (!isTRUE(all(okg %in% names(geneRanges))))  {
            warning("geneRanges does not include names/ranges for all colnames of mgr fflist")
            okg = intersect( okg, names(geneRanges))
            }
        if (!isTRUE(all(oks %in% srids)))  {
            warning("snpRanges does not include names/ranges for all rownames of mgr fflist")
 # match and numeric indexing much more efficient than names for large objects
            oks = match(oks, srids, nomatch=0)
            if (length(oks) == 0) stop("snpRanges does not include any snps in mgr")
            }
 # can't use zero indices in subsetting GRanges...
        ol = findOverlaps(snpRanges[oks[oks>0]], geneRanges[okg] + 
            radius)
        matm = matchMatrix(ol)
        if (nrow(matm) > 0) {
            mgr@fflist[[1]][matm] = 0
        }
    }

transScores = function (smpack, snpchr = "chr1", rhs, K = 20, targdirpref = "tsco", 
    geneApply = lapply, chrnames = paste("chr", as.character(1:22), sep=""), 
    geneRanges = NULL, snpRanges = NULL, radius = 2e+06, renameChrs=NULL, 
    probesToKeep=NULL, batchsize=200, genegran=50, shortfac=10, wrapperEndo=NULL)
{
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
    if (length(chrnames) < 2) 
        stop("must have length(chrnames) >= 2")
    theCall = match.call()
#
# get an image of the expression+genotype data for SNP on specific chromosome snpchr
#
    sms = getSS(smpack, snpchr, renameChrs=renameChrs, probesToKeep=probesToKeep, wrapperEndo=wrapperEndo)
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
    inimgr = eqtlTests(sms[probeId(pnameList[[chrnames[1]]]),   # start the sifting through transcriptome
        ], rhs, targdir = targdir, runname = paste("tsc_", chrnames[1],  # testing on genes in chrom 1
        sep = ""), geneApply = geneApply, saveSummaries = FALSE, genegran=genegran, shortfac=shortfac)
    if (snpchr == chrnames[1]) {
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
        nxtmgr = eqtlTests(sms[probeId(pnameList[[chrnames[j]]]), 
            ], rhs, targdir = targdir, runname = paste("tsctmp", 
            j, sep = ""), geneApply = geneApply, saveSummaries = FALSE,
            genegran=genegran, shortfac=shortfac)
        if (snpchr == chrnames[j]) {
            if (is.null(geneRanges) || is.null(snpRanges)) 
                stop("ranges must be supplied to exclude cis tests")
            cisZero(nxtmgr, snpRanges, geneRanges, radius)
            }
        nxtKinds = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "indscratch.ff", sep = ""), feat = "geneind", ginds = genemap[[j]], 		batchsize=batchsize)
        nxtKscores = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
            "scoscratch.ff", sep = ""), feat = "score", ginds = genemap[[j]],
                batchsize=batchsize)
        updateKfeats(topKscores, nxtKscores, topKinds, nxtKinds, batchsize=batchsize)  
        unlink(filename(nxtmgr@fflist[[1]]))   # kill off scratch materials
        unlink(paste(targdir, "indscratch.ff", sep = ""))
        unlink(paste(targdir, "scoscratch.ff", sep = ""))
    }
    baseout = list(scores = topKscores, inds = topKinds, guniv = guniv, K=K,
	smsanno = annotation(sms),
        snpnames = rownames(inimgr@fflist[[1]]), call = theCall, date=date(), shortfac=shortfac)
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
      scos = cbind(rsco1[i1:i2,], rsco2[i1:i2,])
      ginds = cbind(rind1[i1:i2,], rind2[i1:i2,])
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


bindSnpRanges2mgr = function (mgr, snpRanges, badstart=-2, badend=-1)
{
#
# this will get a GRanges instance from snpRanges congruent to mgr snp info
#
    dimnt = dimnames(mgr@fflist[[1]])
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
    IRanges:::end(fullsr[which(matup > 0)]) = IRanges:::end(snpRanges[matup[matup>0]])
    IRanges:::start(fullsr[which(matup > 0)]) = IRanges:::start(snpRanges[matup[matup>0]])
    fullsr
}

bindGeneRanges2mgr = function (mgr, geneRanges, badstart=-6, badend=-5)
{
#
# this will get a GRanges instance from geneRanges congruent to mgr gene info
#
    dimnt = dimnames(mgr@fflist[[1]])
    genesInMgr = dimnt[[2]]
    fullgr = IRanges(rep(badstart, length(genesInMgr)), rep(badend, length(genesInMgr)))
    sn = seqnames(geneRanges)[1]
    fullgr = GRanges(seqnames = sn, fullgr)
    names(fullgr) = genesInMgr
    matup = match(names(fullgr), names(geneRanges), nomatch = 0)
    IRanges:::end(fullgr[which(matup > 0)]) = IRanges:::end(geneRanges[matup[matup>0]])
    IRanges:::start(fullgr[which(matup > 0)]) = IRanges:::start(geneRanges[matup[matup>0]])
    fullgr
}


cisZero = function (mgr, snpRanges, geneRanges, radius) 
{
#
# use geneRanges+radius to define cis to gene intervals
# find tests of SNPs in these intervals
# set their scores in mgr to zero
#
    realSnpRanges = bindSnpRanges2mgr( mgr, snpRanges )
    realGeneRanges = bindGeneRanges2mgr( mgr, geneRanges )
    ol = findOverlaps(realSnpRanges, realGeneRanges +
        radius)
    matm = matchMatrix(ol)
    if (nrow(matm) > 0) {
        mgr@fflist[[1]][matm] = 0
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



transTab = function( x, snpchr ) {
 gannopkname = x$smsanno
 K = x$K
 sids = rep(x$snpnames, each=K )
 thescos = as.numeric(t(as.ram(x$sco)))/x$shortfac
 theinds = as.numeric(t(as.ram(x$inds)))
 require( gannopkname, character.only=TRUE)
 ganno = gsub(".db", "", gannopkname )
 gchr = sapply(mget(x$guniv, get(paste(ganno, "CHR", sep="")), ifnotfound=NA), "[", 1)
 gsym = sapply(mget(x$guniv, get(paste(ganno, "SYMBOL", sep="")), ifnotfound=NA), "[", 1)
 gent = sapply(mget(x$guniv, get(paste(ganno, "ENTREZID", sep="")), ifnotfound=NA), "[", 1)
 gn = x$guniv[ theinds ]
 gchr = gchr[ theinds ]
 gsym = gsym[ theinds ]
 gent = gent[ theinds ]
 data.frame(snp=sids, sumchisq=thescos, probeid=gn , probechr=gchr,
    sym=gsym, entrez=gent)
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
