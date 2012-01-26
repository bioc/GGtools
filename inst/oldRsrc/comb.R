
setClass("eqtlFDRSummary", representation(
  allpermtops="numeric", obsrd="GRanges", calls="list", nperm="numeric",
  theCall="call", sess="ANY", genome="character", cisRadius="numeric",
  nat05="numeric", nat01="numeric", gene2snpList="list", exTransform="function"))

setGeneric("fullreport", function(x) standardGeneric("fullreport"))
setMethod("fullreport", "eqtlFDRSummary", function(x)
 x@obsrd)

setMethod("show", "eqtlFDRSummary", function(object) {
 nr = length(object@obsrd)
 cat("eqtlTestsSummary instance.  Best observed features (of ", nr, 
    "):\n", sep="")
 todo = min(c(nr, 4))
 allans = object@obsrd
 toshow = object@obsrd[1:todo]
 print(data.frame(probe=names(toshow), 
  chr=as.character(seqnames(toshow)),
  codeStart=start(toshow),
  codeEnd=end(toshow), rsid=elementMetadata(toshow)$rsid,
  snploc=elementMetadata(toshow)$snploc, 
  `chisq(1)`=elementMetadata(toshow)$chisq, check.names=FALSE))
#  `-log10fdr`=elementMetadata(toshow)$score, check.names=FALSE))
 cat("Plug-in FDR based on", object@nperm, "permutations.\n")
 nat01 = sum(elementMetadata(allans)$score >= 2)
 nat05 = sum(elementMetadata(allans)$score >= -log10(0.05))
 cat("There were", nat01, "calls at FDR=0.01;", nat05, "at FDR=0.05.\n", sep=" ")
 cat("call was:\n")
 print(object@theCall)
 cat("===\n")
 cat("use fullreport() for information on all probes tested.\n")
})

setMethod("c", "eqtlFDRtab", function(x, ..., recursive=FALSE) {
 if (recursive) stop("recursive mode not implemented")
 args = unname(list(x, ...))
 cls = sapply(args, class)
 if (!all(cls[1] == cls)) stop("all args to c for eqtlFDRtab must have same class")
 alp = unlist(lapply(args, function(x) x$sorted.all.permtops))
 rd = do.call(c, lapply(args, function(x) x$gro))
 np = do.call(c, lapply(args, function(x) x$nperm))
 g2sls = do.call(c, lapply(args, function(x)x$gene2snpList))
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
 nat01 = sum(elementMetadata(rd)$score >= 2)
 nat05 = sum(elementMetadata(rd)$score >= -log10(0.05))
 new("eqtlFDRSummary", allpermtops=alp, obsrd=rd, calls=calls, nperm=nperm,
          nat01=nat01, nat05=nat05, gene2snpList=g2sls)
})

setMethod("c", "eqtlFDRSummary", function(x, ..., recursive=FALSE) {
 if (recursive) stop("recursive mode not implemented")
 thecall = match.call()
 args = unname(list(x, ...))
 cls = sapply(args, class)
 if (!all(cls[1] == cls)) stop("all args to c for eqtlFDRSummary must have same class")
 np = do.call(c, lapply(args, function(x) x@nperm))
 if (var(np)>0) stop("nperm must be constant over all summaries")
 nperm = np[1]
 alp = unlist(lapply(args, function(x) x@allpermtops))
 rd = do.call(c, lapply(args, function(x) x@obsrd))
 rd = rd[order(elementMetadata(rd)$chisq, decreasing=TRUE)]
 obstops = elementMetadata(rd)$chisq
 sfdr = sapply(obstops, function(x) (sum(alp > x)/nperm)/max(c(1,sum(obstops>x))))
 calls = do.call(c, lapply(args, function(x) x@calls))
 g2sls = do.call(c, lapply(args, function(x) x@gene2snpList))
 lsfdr = ifelse(sfdr==0, 16, -log10(sfdr))
 elementMetadata(rd)$score = lsfdr
 nat01 = sum(elementMetadata(rd)$score >= 2)
 nat05 = sum(elementMetadata(rd)$score >= -log10(0.05))
 new("eqtlFDRSummary", allpermtops=alp, obsrd=rd, calls=calls, nperm=nperm,
          nat01=nat01, nat05=nat05, theCall=thecall, gene2snpList=g2sls)
})
