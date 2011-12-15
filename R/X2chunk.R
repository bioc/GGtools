X2chunk = function( mgr, ffind, start, end, snplocs, anno, useSym=TRUE ) {
  if (!is(ffind, "numeric")) stop("ffind must be numeric")
  if (length(start)>1) stop("start must be atomic")
  lim = IRanges(start,end)
  if (is(mgr, "multffManager")) {
           FFL = mgr$fflist
           SF = mgr$shortfac
           }
  else if (is(mgr, "eqtlTestsManager")) {
           FFL = fflist(mgr)
           SF = shortfac(mgr)
           }
  else stop("mgr must inherit from multffManager or eqtlTestsManager")
  snpNames = rownames(FFL[[ffind]])
  gNames = colnames(FFL[[ffind]])
    #print(ranges(snplocs)[[1]])
    #print(lim)
    indsWithin = function(gr, lim) {
      which(start(gr)>=start(lim) & end(gr)<= end(lim))
    }
    #print(summary(indsWithin(ranges(snplocs)[[1]], lim)))
    #print(summary(ranges(snplocs)[[1]] %in% lim))
    #sn2keep_inds = which(IRanges::"%in%"(ranges(snplocs)[[1]], lim))
    sn2keep_inds = indsWithin( ranges(snplocs)[[1]], lim)
  if (length(sn2keep_inds)==0) stop("start:end too restrictive, no SNP")
  snpInRange = snplocs[ sn2keep_inds, ]
  snp2get = snpInRange[ order(start(snpInRange)), ]$name
  snp2get = intersect(snp2get, snpNames)
  gr = geneRanges( gNames, anno )
    #g2keep_inds = which(IRanges::"%in%"(ranges(gr)[[1]], lim))
    g2keep_inds = indsWithin(ranges(gr)[[1]], lim)
  g2get = gr[g2keep_inds, ]$name
  ans = as.ram(FFL[[ffind]][snp2get, g2get, drop=FALSE])/SF
  if (useSym) {
    gs = geneSyms(colnames(ans), anno)
    if (any(is.na(gs))) gs[which(is.na(gs))] = colnames(ans)[which(is.na(gs))]
    colnames(ans) = gs
  }
  ans
}
  
