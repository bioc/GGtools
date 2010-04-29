
X2chunk = function( mgr, ffind, start, end, snplocs, anno, useSym=TRUE ) {
  if (!is(ffind, "numeric")) stop("ffind must be numeric")
  if (length(start)>1) stop("start must be atomic")
  lim = IRanges(start,end)
  snpNames = rownames(mgr$fflist[[ffind]])
  gNames = colnames(mgr$fflist[[ffind]])
    sn2keep_inds = IRanges::"%in%"(ranges(snplocs)[[1]], lim)
  if (length(sn2keep_inds)==0) stop("start:end too restrictive, no SNP")
  snpInRange = snplocs[ sn2keep_inds, ]
  snp2get = snpInRange[ order(start(snpInRange)), ]$name
  snp2get = intersect(snp2get, snpNames)
  gr = geneRanges( gNames, anno )
    g2keep_inds = IRanges::"%in%"(ranges(gr)[[1]], lim)
  g2get = gr[g2keep_inds, ]$name
  ans = as.ram(mgr$fflist[[ffind]][snp2get, g2get])/mgr$shortfac
  if (useSym) {
    gs = geneSyms(colnames(ans), anno)
    if (any(is.na(gs))) gs[which(is.na(gs))] = colnames(ans)[which(is.na(gs))]
    colnames(ans) = gs
  }
  ans
}
  

