
setClass("eqtlFDRSummary", representation(
  allpermtops="numeric", obsrd="GRanges", calls="list",
  theCall="call", sess="ANY", genome="character", cisRadius="numeric",
  nperm="numeric"))

setMethod("show", "eqtlFDRSummary", function(object) {
 nr = length(object@obsrd)
 cat("eqtlTestsSummary instance.  Best observed features (of", nr, 
    "):\n", sep="")
 todo = min(c(nr, 3))
 print(object@obsrd[1:todo])
})

setMethod("c", "eqtlFDRtab", function(x, ..., recursive=FALSE) {
 if (recursive) stop("recursive mode not implemented")
 args = unname(list(x, ...))
 cls = sapply(args, class)
 if (!all(cls[1] == cls)) stop("all args to c for eqtlTestsManager must have same class")
 alp = unlist(lapply(args, function(x) x$sorted.all.permtops))
 rd = do.call(c, lapply(args, function(x) x$gro))
 np = do.call(c, lapply(args, function(x) x$nperm))
 if (!(all(np[1] == np))) warning("different numbers of permutations in different eqtlFDRtab, using max")
 nperm = max(np)
 calls = lapply(args, function(x) x$thecall)
 rd = rd[ order(elementMetadata(rd)$score, elementMetadata(rd)$chisq, decreasing=TRUE), ]
 obstops = elementMetadata(rd)$chisq
 sfdr = sapply(obstops, function(x) (sum(alp > x)/nperm)/max(c(1,sum(obstops>x))))
 nn = names(elementMetadata(rd))
 nn[nn=="score"] = "chrWiseScore"
 names(elementMetadata(rd)) = nn
 lsfdr = ifelse(sfdr==0, 16, -log10(sfdr))
 elementMetadata(rd)$score = lsfdr
 new("eqtlFDRSummary", allpermtops=alp, obsrd=rd, calls=calls)
})
