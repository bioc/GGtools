best.trans.eQTLs = function(smpack, rhs, genechrnum, snpchrnum,
	K = 20, targdirpref="tsco", batchsize=200, radius=2e6,
        genequeryprefix="", snploadprefix="chr", snplocprefix="chr", geneannopk,
	snpannopk, exFilter=function(x)x, smFilter=function(x)x,
	geneApply=lapply) {
   transScores( smpack=smpack, snpchr=paste( snploadprefix, snpchrnum, sep=""), 
		rhs=rhs, K=K, 
		targdirpref=targdirpref, chrnames=genechrnum,
		gchrpref=genequeryprefix, schrpref=snplocprefix,
		radius=radius, shortfac=10, wrapperEndo=smFilter,
		geneannopk=geneannopk, snpannopk=snpannopk )
}
		
#meta.best.trans.eQTLs = function(smpackvec, rhslist, genechrnum, snpchrnum,
#	K = 20, targdirpref="multtsco", batchsize=200, radius=2e6,
#        genequeryprefix="", snploadprefix="chr", snplocprefix="chr", geneannopk,
#	snpannopk, exFilterList=function(x)x, smFilterList=function(x)x,
#	geneApply=lapply) {
#   mtransScores( smpackvec, rhslist=rhslist, chrnames=genechrnum, snpchr=snpchrnum,
#         targdirpref=targdirpref, geneApply=geneApply, geneannopk=geneannopk,
#         snpannopk=snpannopk, radius=radius,  
#
#mtransScores = function (smpackvec, snpchr = "chr1", rhslist, K = 20, targdirpref = "multtsco", 
#    geneApply = lapply, chrnames = paste("chr", as.character(1:22), sep=""), 
#    geneRanges = NULL, snpRanges = NULL, radius = 2e+06, renameChrs=NULL,
#    batchsize=200, genegran=50, probesToKeep=NULL, shortfac=10, wrapperEndo=NULL,
#    gchrpref=genequeryprefix, ) 
#{
##
## objective is a small-footprint accumulation of trans-eQTL tests 
##      summed over different cohorts
##  smpackvec is a vector of lightweight smlSet package name
##  snpchr is the chromosome for which SNPs will be tested
##  rhslist is the right hand side of formula for snp.rhs.tests in snpTests
##  K is the number of best features to be retained as we explore the transcriptome
##
##
#    if (length(chrnames) < 2) 
#        stop("must have length(chrnames) >= 2")
#    theCall = match.call()
#    require(GGtools)
##
## get an image of the expression+genotype data for SNP on specific chromosome snpchr
##
#    smsl = lapply(smpackvec, function(x) getSS(x, paste(snploadprefix, snpchr, sep="")))
#    smsl = lapply(1:length(smsl), function(z) exFilterList(z))  # essential for dealing with multiple cell types
#    smsl = lapply(1:length(smsl), function(z) smFilterList(z))
#
#    smsl = makeCommonSNPs(smsl) # could be optional
#    names(smsl) = smpackvec
#    sms = smsl[[1]]
#    guniv = featureNames(sms)   # universe of probes
#    smanno = gsub(".db", "", annotation(sms))
#    require(paste(smanno, ".db", sep = ""), character.only = TRUE)
#    clcnames = gsub("chr", "", chrnames)  # typical chrom nomenclature of bioconductor
#    pnameList = mget(clcnames, revmap(get(paste(smanno, "CHR", 
#        sep = ""))))
# # be sure to use only genes that are on arrays in sms
#    pnameList = lapply(pnameList, function(x) intersect(x, guniv))
#    names(pnameList) = chrnames  # now the universe of probes is split into chromsomes
#    todrop = which(sapply(pnameList, length)==0)
#    if (length(todrop)>0) pnameList = pnameList[-todrop]
#    genemap = lapply(pnameList, function(x) match(x, guniv))  # numerical indices for probes
#    nchr_genes = length(names(pnameList))
#    targdir = paste(targdirpref, snpchr, sep="")
#    cursmsl = lapply(smsl, function(x) x[probeId(pnameList[[chrnames[1]]]),])
#    inimgr = meqtlTests(cursmsl,   # start the sifting through transcriptome
#         rhslist, targdir = targdir, runname = paste("tsc_", chrnames[1],  # testing on genes in chrom 1
#        sep = ""), geneApply = geneApply, shortfac=shortfac)
#    rm(cursmsl); gc()
# # build the map for current snpchr
#    mapobj = getCisMap( radius = radius, gchr = paste(gchrpref, chrnames[j], sep=""),
#                  schr = paste(schrpref, snpchr, sep=""), geneannopk=geneannopk, snpannopk = snpannopk )
## from transScores
#        if (snpchr == chrnames[1]) {
#            cisZero(inimgr, mapobj@snplocs, mapobj@generanges, radius=0)   # if SNP are on chrom 1, exclude cis
#                             # the gene ranges supplied are already augmented by radius
#        }
#
#    topKinds = topKfeats(inimgr, K = K, fn = paste(targdir, "/",  # sort and save
#        snpchr, "_tsinds1_1.ff", sep = ""), feat = "geneind", 
#        ginds = genemap[[1]], batchsize=batchsize)
#    topKscores = topKfeats(inimgr, K = K, fn = paste(targdir, 
#        "/", snpchr, "_tssco1_1.ff", sep = ""), feat = "score", 
#        ginds = genemap[[1]], batchsize=batchsize)
#    unlink(filename(inimgr@fflist[[1]]))
#    for (j in 2:nchr_genes) {    # continue sifting through transcriptome
#        cat(j)
#        gc()
#        cursmsl = lapply(smsl, function(x) x[probeId(pnameList[[chrnames[j]]]),])
#        nxtmgr = meqtlTests(cursmsl, rhslist, targdir = targdir, 
#            runname = paste("tsctmp", 
#            j, sep = ""), geneApply = geneApply, shortfac=shortfac)
#        rm(cursmsl); gc()
#        if (snpchr == chrnames[j]) {
#            if (is.null(geneRanges) || is.null(snpRanges)) 
#                stop("ranges must be supplied to exclude cis tests")
#            cisZero(nxtmgr, snpRanges, geneRanges, radius)
#            }
#        nxtKinds = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
#            "indscratch.ff", sep = ""), feat = "geneind", ginds = genemap[[j]], batchsize=batchsize)
#        nxtKscores = topKfeats(nxtmgr, K = K, fn = paste(targdir, 
#            "scoscratch.ff", sep = ""), feat = "score", ginds = genemap[[j]], batchsize=batchsize)
#        updateKfeats(topKscores, nxtKscores, topKinds, nxtKinds, batchsize=batchsize)  
#        unlink(filename(nxtmgr@fflist[[1]]))   # kill off scratch materials
#        unlink(paste(targdir, "indscratch.ff", sep = ""))
#        unlink(paste(targdir, "scoscratch.ff", sep = ""))
#    }
#
#    baseout = list(scores = topKscores, inds = topKinds, guniv = guniv, K=K, snpchr=snpchr,
#        chrnames=chrnames,
#	smsanno = annotation(sms),
#        snpnames = rownames(inimgr@fffile), call = theCall, date=date(), shortfac=shortfac)
#    new("transManager", base=baseout)
#}
