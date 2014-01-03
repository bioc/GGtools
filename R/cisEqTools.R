
ciseqByCluster = function( cl, pack = "yri1kgv",
      outprefix = "yrirun",
      finaltag = "partyri100k",
      chromsToRun = 1:22,  # if length is C will use 2C nodes 
      targetfolder = "/freshdata/YRI_3",
      radius = 100000L,
      nperm = 3L,
      ncoresPerNode = 8,
      numPCtoFilter = 10,  # filtering actions can be made more generic, 1/1/2014
      lowerMAF = .02,
      geneannopk = "lumiHumanAll.db",
      snpannopk = "SNPlocs.Hsapiens.dbSNP.20120608",
      smchrpref = "chr",
      tmpForSort = "/tmp",
      numtiles = 200,
      postProcCores=12) {
#
# this program will split all chromsomes in roughly equal halves (by probe) and
# compute All.cis for all probes ... snps may recur in the different halves
#
# performance: assumes parallel-based cluster instance defines major nodes, each is 
# treated as multicore.  this should be runnable on a scalar machine by setting cl
# to be a single PSOCK cluster at "localhost", and setting ncoresPerNode = 1
#
      #nodeNames = rep("localhost", numNodes),
      #numNodes = 8,
#
  stopifnot(inherits(cl, "cluster")) # can be any sort of cluster instance from parallel package
  numNodes = length(cl)  # assumes is.list(cl) is true ...
  chkNcores = unlist(clusterApply(cl, 1:numNodes, function(x) {library(parallel);  detectCores()}))
  if (any(chkNcores < ncoresPerNode)) warning("some nodes have fewer than ncoresPerNode cores")
  if (!file.exists(targetfolder)) try(system(paste0("mkdir ", targetfolder)))
  if (!file.exists(targetfolder)) stop(paste0("cannot  create ", targetfolder))
  library(parallel)
  firstHalf <<- function(x) x[1:floor(length(x))/2]
  secondHalf <<- function(x) x[-(1:floor(length(x))/2)]
  setupSplit = function(nodeset=1:numNodes) {
     clusterApply(cl, nodeset, function(w) {
      library(parallel)  # get resources
      library(GGtools)
      library(pack, character.only=TRUE)
      library(geneannopk, character.only=TRUE)
      library(snpannopk, character.only=TRUE)
      #library(Rsamtools)  # need for export but CMD check complains
      #library(rtracklayer)
      p1 = "Rsamtools"
      p2 = "rtracklayer"
      library(p1, character.only=TRUE)
      library(p2, character.only=TRUE)
      cc = new("CisConfig")  # configure search, except for choice of chromosome
      smpack(cc) = pack
      nperm(cc) = nperm
      geneannopk(cc) = geneannopk
      radius(cc) = radius
      smchrpref(cc) = smchrpref
      geneApply(cc) = mclapply  # so genes are dispatched to cores
      options(mc.cores=ncoresPerNode)
      cc <<- cc
     } )
  }
  setupSplit(1:numNodes)  # causes library attachment on all nodes and generation of CisConfig instances there
  
  runOneSplit <<- function(chrtag, suffix="A", splitter=firstHalf, gffOnly=TRUE) { # assumes cc is defined locally as the config
    if (!exists("cc")) stop("did not find cc for local CisConfig manipulation")
    chrnames(cc) = as.character(chrtag)
    folderStem(cc) = paste0(folderStem(cc), "_", chrtag, suffix)
    smFilter(cc) = function(x) {   # late filtering
        fn = featureNames(x);
        if (numPCtoFilter > 0) tmp = MAFfilter( clipPCs(x, 1:numPCtoFilter), lower=lowerMAF )
        else tmp = MAFfilter( x, lower=lowerMAF )
        tmp[probeId(splitter(fn)),]
        }
    obn = paste0(outprefix, "_", chrnames(cc), suffix)
    fn = paste0(obn, ".rda")
    res <- try(All.cis(cc))
    if (inherits(res, "try-error")) return(res)

          getmindist = function(snplocs, probes, geneannopk="lumiHumanAll.db") {
           require(geneannopk, character.only=TRUE, quietly=TRUE)
           stub = gsub(".db$", "", geneannopk)
           locenv = get(paste0(stub, "CHRLOC"))
           locendenv = get(paste0(stub, "CHRLOCEND"))
           gloc = abs(sapply(mget(probes, locenv), "[", 1))
           gend = abs(sapply(mget(probes, locendenv), "[", 1))
           ifelse(
                    snplocs <= gend & snplocs >= gloc, 0,
                    ifelse( snplocs > gend, snplocs-gend, gloc-snplocs ) )
          }
    res$mindist = getmindist(res$snplocs, res$probeid, geneannopk=geneannopk)
    assign(obn, res)

    if(!gffOnly) {
       save(list=obn, file=fn)    # save to local disk
       system(paste0("cp ", fn, " ", targetfolder ))  # copy to target folder
       }
#
#
# define transformer
  cr2gff = function(cr, obn, targetfolder=".") {
      p1 = "GenomicRanges"
      p2 = "Rsamtools"
      require(p1, character.only=TRUE)
      require(p2, character.only=TRUE)
      res = as(cr, "GRanges")
      sl = IRanges(res$snplocs, width=1)
      ranges(res) = sl
      names(res) = NULL
      sn = seqnames(res)
      sn = gsub("chr", "", sn)
      o = order( as.numeric(sn), start(res) )
      res = res[o]
      gffFile = paste0(targetfolder, "/", obn, ".gff3")
  #    seqlevels(res) = gsub("chr", "", seqlevels(res))
      export.gff3( res, gffFile )
 #     owd = getwd()
 #     setwd(targetfolder)
 #     bgzip( gffFile , overwrite=TRUE)
 #     indexTabix( paste0(gffFile, ".gz") , format="gff")
 #     setwd( owd )
  }
#
# continue processing to tabix-indexed gff3 based on SNP address
#
    cr2gff( as(res, "GRanges"), obn, targetfolder )
    NULL
  }
  clusterExport(cl, "runOneSplit")
  clusterExport(cl, "firstHalf")
  clusterExport(cl, "secondHalf")
 
  njobs = 2*length(chromsToRun)
  chrtags = as.character(rep(chromsToRun, each=2))
  
  reqlist = vector("list", njobs)
  j <- 1
  for (i in 1:length(chromsToRun)) {
   reqlist[[j]] = list( chr=chrtags[j], tag="A", splitter=firstHalf )  # need to distinguish splitter elements
   reqlist[[j+1]] = list( chr=chrtags[j+1], tag="B", splitter=secondHalf )  # therefore loop is complex
   j <- j+2
  }
  reqlist <<- reqlist
  
  clusterApplyLB(cl, reqlist, function(x) runOneSplit(x[["chr"]], x[["tag"]], x[["splitter"]]))

curd = getwd()
on.exit(setwd(curd))
setwd(targetfolder)
gffprocess(finaltag, n_in=njobs, headpatt=paste0("_", chrtags[1], "A"), tmpForSort=tmpForSort)
myti = simpleTiling(numtiles)
pmti = myti[seqnames(myti) %in% paste0("chr", chromsToRun)]
require(foreach)
cgff2dt(paste0(finaltag, ".gff3.gz"), pmti)

}
