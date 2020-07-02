
suppressPackageStartupMessages(library(GGtools))

# configure
   cc = new("CisConfig")
   chrnames(cc) = "21"
   genome(cc) = "hg19"
   nperm(cc) = 2L
   lkp = try(library(parallel))
   if (!inherits(lkp, "try-error")) {
      nc = min(10, detectCores())
      options(mc.cores=nc)
      geneApply(cc) = mclapply
      }
   estimates(cc) = FALSE
   set.seed(1234)
#   system.time(f1 <- cisScores( cc ))
 #
 # demonstrate adding annotation on chromatin state and gwas status
 #
 eprops = function(ans) {
 #
 # only adds fields to values() of the input
 #
  data(hmm878)
  ac = as.character
  eqr = GRanges(ac(seqnames(ans)), IRanges(ans$snplocs, width=1))
  fo = findOverlaps(eqr, hmm878)
  chromcat878 = factor(rep("none", length(ans)), levels=c(unique(hmm878$name), "none"))
  chromcat878[ queryHits(fo) ] = factor(hmm878$name[subjectHits(fo)])
  ans$chromcat878 = chromcat878
 
  if (require(gwascat)) {
    data(gwastagger)
    isgwashit = 1*(overlapsAny(eqr, gwastagger) | ans$snp %in% gwastagger$tagid) # allow match by loc or name
    ans$isgwashit = isgwashit
    }
  ans
 }
 extraProps(cc) = eprops
 set.seed(1234)
 rhs(cc) = ~1-1
if (.Platform$OS.type != "windows") {
 (f2 <- cisScores( cc ))
 isTRUE(sum(f2$fdr < 0.05) > 0)  # can change with annotation or location changes, check serialized results if necessary
}
TRUE
