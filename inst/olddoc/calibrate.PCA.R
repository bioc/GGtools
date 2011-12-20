calibrate.PCA = function(smls, indices=1:200, npcs=c(1,2,5,10,20,40), ffind=1, 
   gapplier=lapply) {
#
# program computes PCA on all genes/samples, then restricts testing to genes in 'indices'
#
# parameters
#   smls = smlSet -- should represent whole array as basis for PCA
#   indices = some indexing into features to reduce the testing problem -- could
#      limit to one chromosome or more sharply for example
#   npcs = number of PC to put as independent variables in candidate models
#   ffind = which fflist element is to be used to get p-values
#   gapplier = how to iterate over genes, could be mclapply
#
  if (max(npcs) > length(featureNames(smls))) stop("not enough genes to compute requested number of PCs")
  if (length(smList(smls))>1) stop("get a single chrom smlset for now")
  pc = prcomp(t(exprs(smls)))$x
  rownames(pc) = sampleNames(smls)
  pData(smls) = cbind(pData(smls), pc)
  fmlas = lapply(npcs, function(x) 1:x)
  fmlas = lapply(fmlas, function(x) paste("PC", x, sep=""))
  fmlas = lapply(fmlas, paste, collapse="+")
  fmlas = lapply(fmlas, function(x) paste("~", x, sep=""))
  fmlas = lapply(fmlas , as.formula)
  ans = list()
  smls = smls[indices,]
  allprobes = featureNames(smls)
  for (i in 1:length(fmlas)) {
    try(system("rm -rf foo"))
    t1 = eqtlTests( smls, fmlas[[i]], geneApply=gapplier )
    ans[[i]] = sapply(allprobes, function(x) topFeats(probeId(x), mgr=t1, ffind=ffind, n=1))
  }
  ans = lapply(ans, function(x) 1-pchisq(x, t1@df))
  names(ans) = paste("#PC=", npcs, sep="")
  ans
}

#
# we demonstrate with test data inside GGtools package
#
library(GGtools)
data(hmceuB36.2021)
hlit = hmceuB36.2021[ chrnum("20"), ]
library(illuminaHumanv1.db)
g20 = get("20", revmap(illuminaHumanv1CHR))
g20 = intersect(g20, featureNames(hmceuB36.2021))
library(multicore)
options(cores=12)
dd = calibrate.PCA( hlit, probeId(g20[1:150]), gapplier=mclapply )
#
# following code makes the kind of plot seen before
#
MM = sapply(dd, function(x) sapply(c(1e-3,1e-4,1e-5,1e-6,1e-7), function(z)sum(x<z))) 
matplot(c(1,2,5,10,20,40), t(MM), type="l", xlab="#PC", ylab="num probes w/p< legend")
legend(1,140, col=1:5, lty=1:5, legend=c(paste("10^-", 3:7, sep="")))

