
#> fullceu100k_dtplus[1,]
#   seqnames start   end width strand      source             type score phase
#1:     chr1 13302 13302     1      + rtracklayer sequence_feature  1.99    NA
#           snp snplocs  ests   se       fdr            probeid       MAF
#1: rs180734498   13302 -0.02 0.01 0.8185385 B7ukGTu2k86e0Tqg_k 0.1835443
#   dist.mid permScore_1 permScore_2 permScore_3 mindist    oldfdr chromcat878
#1: -56247.5        0.41        0.44        0.36   55789 0.8223583 11_Weak_Txn
#   isgwashit       distfac    fdrfac     maffac
#1:         0 (5e+04,1e+05] (0.2,1.1] (0.1,0.25]
#

update_fdr_filt = function(tab, 
    filt=function(x)x, by=c("pairs", "snps", "probes")[1]) {
#
# take a table, a filtering function, and a scope specification
# for FDR -- 'pairs' implies every probe-snp pair is a test
# 'snps' implies we look at all probes cis to snp and use
#   strongest association as score for the snp
# 'probes' implies we look at all snps cis to the probe and use
#   strongest association as score for the probe
#
 require(GGtools, quietly=TRUE)
 tab = filt(tab)
 psinds = grep("permScore", colnames(tab), value=TRUE)
 #ptab = as.data.frame(tab[,psinds,with=FALSE]) # decent
 #pscores = as.numeric(unlist(ptab))  # slow
 nr = nrow(tab)
 pscores = vector("numeric", nr*length(psinds))
 for (np in 1:length(psinds)) 
   pscores[ (((np-1)*nr)+1):(np*nr) ] = tab[[psinds[np]]]
 if (by == "pairs") {
   newfd = pifdr( tab$score, pscores )
   }
 else {
   if (by == "snps") byvar = "snp"
   else if (by == "probes") byvar = "probeid"
   base = tab[, max(score), by=byvar]
   maxbysnp = base$V1
   ol = list()
   pnames = grep("permScore", names(tab))  # numeric indices
   for (i in 1:length(pnames))  
     {
     tab$score = tab[,pnames[i],with=FALSE]   # force reuse of score field
     ol[[i]] = tab[, max(score), by=byvar]$V1
     }
   newfd = pifdr( maxbysnp, as.numeric(unlist(ol)))
   tab = base
   }
 tab$fdr = newfd
 tab
}
  
# sensitivity analysis for eQTL search on MAF and distance of search

filtgen.maf.dist = function(maf.dist, 
    validate.tab=function(tab)all(c("mindist", "MAF") %in% colnames(tab))) {
  stopifnot(is.atomic(maf.dist))
  stopifnot(length(maf.dist)==2)
  maf = maf.dist[1]
  dist = maf.dist[2]
#
# illustrative closure generator: take a data.table and filter on two predicates
# take rows with MAF >= maf and mindist <= dist
#
  function(tab) {
     stopifnot(isTRUE(validate.tab(tab)))
     tab[ tab$mindist <= dist & tab$MAF >= maf, ]
     }
  }

eqsens_dt = function(dt, filtgen=filtgen.maf.dist, 
   by=c("pairs", "snps", "probes")[1],
   targfdrs=c(.05, .01, .005), parmslist=list(
   mafs = c( .025, .05, .075, .1, .125 ),
   dists = c( 1000, 5000, 10000, 25000, 50000, 100000 ) ) ) {
#
# returns counts
#
  parmset = data.matrix(do.call(expand.grid, parmslist))
  ntune = nrow(parmset)
  ans = foreach( curp=1:ntune ) %dopar% {
    tmp = update_fdr_filt( tab=dt, filt=filtgen( parmset[curp, ] ), by=by )
    sapply(targfdrs, function(x)sum(tmp$fdr <= x))
    }
  hold = t(sapply(ans,force))
  colnames(hold) = paste0("at_", targfdrs)
  cbind(parmset, hold)
}

#pairs.df = data.frame(maf=allarg[,1], dist=allarg[,2], at05=NA, at01=NA, at005=NA )
#
## counting at various FDR
#
#ssm(library(GGtools))
#source("pureSensSoft.R")
#load("partceu100k_dt.rda")
#
#library(parallel)
#options(mc.cores=8)
#
#bag = mclapply( 1:nrow(pairs.df), function(x) {
#  tmp = update_fdr_filt( partceu100k_dt,
#       filtgen.maf.dist(allarg[x,1], allarg[x,2] ), by="probes" )
#  ans = sapply(c(.05, .01, .005), function(x) sum(tmp$fdr <= x) )
#  rm(tmp)
#  gc()
#  ans
#})
#
# 
#
