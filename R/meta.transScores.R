meta.transScores = function (smpackvec = c("GGdata", "hmyriB36"), 
    snpchr = "22", rhsList=list(~1, ~1), K = 20, targdirpref = "mtsco", 
    geneApply = lapply, chrnames = as.character(21:22), 
    radius = 2e+06, renameChrs=NULL, 
    probesToKeep=NULL, batchsize=200, genegran=50, shortfac=10, wrapperEndo=NULL,
    geneannopk = "illuminaHumanv1.db", snpannopk = snplocsDefault(),
    gchrpref = "", schrpref="ch", 
    exFilterList= list(function(x)x, function(x)x),
    SMFilterList = list(function(x)x, function(x)x), SSgen=GGBase::getSS)
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
    #if (length(chrnames) < 2) 
    #    stop("must have length(chrnames) >= 2")
    theCall = match.call()
#
# get an image of the expression+genotype data for SNP on specific chromosome snpchr
#

    cat("get data...")
    smsList = lapply(1:length(smpackvec), function(x) SSgen(
                smpackvec[x], snpchr, exFilter=exFilterList[[x]]))
    cat("run smFilter...")   # run earlier than in eqtlTests
    for (i in 1:length(smsList))
         smsList[[i]] =  SMFilterList[[i]](smsList[[i]])
#   get common featurenames
    cat("get common feature names...")
    allpn = lapply(smsList, featureNames)
    commpn = allpn[[1]]
    for (i in 2:length(allpn)) commpn = intersect(commpn, allpn[[i]])
    smsList = lapply(smsList, function(x) x[probeId(commpn),])  # common probes established

    if (!is.null(renameChrs)) snpchr=renameChrs
    guniv = featureNames(smsList[[1]])   # universe of probes
    annos = sapply(smsList, annotation)
    if (!all(annos == annos[1])) stop("distinct annotation fields for smlSets in smpackvec; should be identical across all")
    smanno = gsub(".db", "", annos[1])
    require(paste(smanno, ".db", sep = ""), character.only = TRUE)
    clcnames = gsub("chr", "", chrnames)  # typical chrom nomenclature of bioconductor
    pnameList = mget(clcnames, revmap(get(paste(smanno, "CHR", 
        sep = ""))))
    plens = sapply(pnameList, length)
    if (all(plens == 0)) stop("chromosome map for probes is empty")
 # be sure to use only genes that are on arrays in sms
    pnameList = lapply(pnameList, function(x) intersect(x, guniv))
    names(pnameList) = chrnames  # now the universe of probes is split into chromsomes
    todrop = which(sapply(pnameList, length)==0)
    if (length(todrop)>0) pnameList = pnameList[-todrop]
    genemap = lapply(pnameList, function(x) match(x, guniv))  # numerical indices for probes
    nchr_genes = length(names(pnameList))
    targdir = paste(targdirpref, snpchr, sep="")
#
#  sanitize smList probes for current chromosome
#
    cursmsList = lapply(smsList, function(x) x[probeId(pnameList[[chrnames[1] ]]) ] )
    
#
# start with first element of chrnames vector
#
#eqtlTests(sms[probeId(pnameList[[chrnames[j]]]), ]...
    inimgr = meqtlTests(listOfSmls=cursmsList,   # start the sifting through transcriptome
        rhslist=rhsList, targdir = targdir, 
        runname = paste("mtsc_", chrnames[1],  sep = ""), 
        geneApply = geneApply, shortfac=shortfac)
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
#  sanitize smList probes for current chromosome
#
       cursmsList = lapply(smsList, function(x) x[probeId(pnameList[[chrnames[j] ]]) ] )
       nxtmgr = meqtlTests(listOfSmls=cursmsList,   # start the sifting through transcriptome
          rhslist=rhsList, targdir = targdir,
          runname = paste("mtsc_", chrnames[j],  sep = ""), 
          geneApply = geneApply, shortfac=shortfac)
        if (gsub("chr", "", snpchr) == gsub("chr", "", chrnames[j])) {
            mapobj = getCisMap( radius = radius, gchr = paste(gchrpref, chrnames[j], sep=""),
                  schr = paste(schrpref, gsub("chr", "", snpchr), sep=""), geneannopk=geneannopk, snpannopk = snpannopk )
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
	smsanno = annotation(smsList[[1]]),
# previous to 21 aug 2012 snpnames were assigned here from rownames inimgr@fffile; these are now
# memoed above
        snpnames = kpsnpnames, call = theCall, date=date(), shortfac=shortfac)
    new("transManager", base=baseout)
}

