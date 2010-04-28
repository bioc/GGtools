
ffchunk = function( mgr, ffind, start, end, snplocs, anno ) {
  if (!is(ffind, "numeric")) stop("ffind must be numeric")
  if (length(start)>1) stop("start must be atomic")
  lim = IRanges(start,end)
  snpNames = rownames(mgr$fflist[[ffind]])
  gNames = colnames(mgr$fflist[[ffind]])
  sn2keep_inds = which( ranges(snplocs)[[1]] %in% lim )
  if (length(sn2keep_inds)==0) stop("start:end too restrictive, no SNP")
  snpInRange = snplocs[ sn2keep_inds, ]
  snp2get = snpInRange[ order(start(snpInRange)), ]$name
  snp2get = intersect(snp2get, snpNames)
  gr = geneRanges( gNames, anno )
  g2keep_inds = which( ranges(gr)[[1]] %in% lim )
  g2get = gr[g2keep_inds, ]$name
  as.ram(mgr$fflist[[ffind]][snp2get, g2get])
}
  

