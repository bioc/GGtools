
geneTrack = function(mgr, gn, chrtag, locdata, dropDups=TRUE, mlog10p=TRUE, minchisq=.001) {
#
# intention is to use ucsc table on snp locations to orient test results for display in
# browser -- locs2keep will have scores and can be exported to wig
#
 snpnames = rownames(mgr$fflist[[chrtag]])
 genenames = colnames(mgr$fflist[[chrtag]])
 if (!(gn %in% genenames)) stop("gene name not found")
 scores = mgr[, probeId(gn)][[chrtag]]
 locnames = locdata$name
 locs2keep = locdata[ which(locnames %in% snpnames), ]
 names(scores) = snpnames
 newsco = scores[ locs2keep$name ]
 locs2keep$score = newsco
 dd = duplicated(start(locs2keep))
 if (any(dd) & dropDups) locs2keep = locs2keep[-which(dd), ]
 if (mlog10p) locs2keep$score = -log10(1-pchisq(locs2keep$score+minchisq, mgr$df))
 locs2keep
}
