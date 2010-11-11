calibrate.PCA = function(smls, npcs=c(1,2,5,10,20,40), ffind=1, pvthresh=c(1/(10^(4:8))),
   gapplier=lapply) {
  if (max(npcs) > length(featureNames(smls))) stop("not enough genes to compute requested number of PCs")
  if (length(smList(smls))>1) stop("get a single chrom smlset for now")
  allprobes = featureNames(smls)
  pc = prcomp(t(exprs(smls)))$x
  rownames(pc) = sampleNames(smls)
  pData(smls) = cbind(pData(smls), pc)
  fmlas = lapply(npcs, function(x) 1:x)
  fmlas = lapply(fmlas, function(x) paste("PC", x, sep=""))
  fmlas = lapply(fmlas, paste, collapse="+")
  fmlas = lapply(fmlas, function(x) paste("~", x, sep=""))
  fmlas = lapply(fmlas , as.formula)
  ans = list()
  for (i in 1:length(fmlas)) {
    try(system("rm -rf foo"))
    t1 = eqtlTests( smls, fmlas[[i]], geneApply=gapplier )
    ans[[i]] = sapply(allprobes, function(x) topFeats(probeId(x), mgr=t1, ffind=ffind, n=1))
  }
  ans
}
