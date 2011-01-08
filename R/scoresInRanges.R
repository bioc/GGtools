upstr = function (gr, rad = 5000)
{
    ee = start(gr) - 1
    ss = start(gr) - rad
    start(gr) = ss
    end(gr) = ee
    gr
}

downstr = function(gr, rad=5000) {
  ss = end(gr)+1
  ee = end(gr)+rad
  start(gr) = ss
  end(gr) = ee
  gr
}
flankingOnly = function(gr, rad=5000) {
  ans = c(upstr(gr,rad), downstr(gr,rad))
  lgr = length(gr)
  reo = as.integer(rbind(1:lgr, (lgr+1):(2*lgr)))
  ans = ans[reo]
  values(ans)$gname = rep(names(gr),each=2)
  ans
}


geneLimits = function( anno="org.Hs.eg.db", chr="20" ) {
 require(anno, character.only=TRUE)
 clnanno = gsub(".db", "", anno)
 gna = get(chr, revmap( get(paste(clnanno, "CHR", sep="")) ) )
 gs = mget(gna, get(paste(clnanno, "CHRLOC", sep="")) ) 
 ge = mget(gna, get(paste(clnanno, "CHRLOCEND", sep="")) ) 
 gs = abs(sapply(gs, "[", 1))
 ge = abs(sapply(ge, "[", 1))
 bad = which(is.na(gs) | is.na(ge))
 ans = GRanges( IRanges(gs[-bad], ge[-bad]), seqnames=paste("chr", chr, sep=""))
 names(ans) = gna[-bad]
 ans
}

extendGR = function( gr, siz=1e6 )
  punion(shift(gr, -siz), shift(gr, siz), fill.gap=TRUE)

containedSnps = function(geneRanges, snpRanges) {
 #
 # returns a list with one element per gene in geneRanges
 # giving the names of the SNP that lie within the gene's range
 #
  mm = findOverlaps( geneRanges, snpRanges )@matchMatrix
  if (nrow(mm) == 0) {
     ans = lapply(1:length(geneRanges), function(x)NA)
     names(ans) = names(geneRanges)
     return(ans)
  }
  sinds = split( mm[,2], mm[,1] )  # split snp indexes by gene index
  ugi = unique(mm[,1])
  gn = names(geneRanges)[ugi]
  ans = lapply( 1:length(sinds), function(x) names(snpRanges)[sinds[[x]] ] )
  names(ans) = gn
  ans
}

#scoresInRanges = function( mgr, geneRanges, snpRanges, applier=lapply ) {
# snm = unlist(lapply(mgr@fflist, rownames))
# iniRangeNames = names(snpRanges)
# onboard = intersect(snm, iniRangeNames)
# mm = match(onboard, iniRangeNames)
# snpRanges = snpRanges[ mm ] # shld be faster than intersect(snm, names(snpRanges)) ]
# snps = containedSnps( geneRanges, snpRanges )
# gn = names(geneRanges)
# applier( 1:length(snps), function(x) mgr[ rsid(snps[[x]]), probeId(gn[x]) ])
#}
 

scoresInRanges = function (mgr, geneRanges, snpRanges, applier = lapply, 
   ffind=NULL, matchProbeNames=TRUE) 
{
 #
 # mgr is an eqtlTestsManager instance
 # geneRanges will typically be a GRanges extended for 'cis'
 # snpRanges can be a general snp metadata GRanges
 #  oct 22-- introduced ffind to reduce scope of mgr to be examined if geneRanges is limited to one chr
 #           also added filtering of geneRanges to genes available in mgr
 # ffind is a detail telling which piece of the mgr (typically chrom) is in use
 # matchProbeNames tells us to restrict attention to ranges in GRanges that
 #   have a name matching a probe in the mgr
 #
# jan 8 2011 -- ffind needs to be set -- could cause problems for 
    if (is.null(ffind)) stop("must set ffind")
    if (!is(geneRanges, "GRanges")) stop("geneRanges must inherit from GRanges")
    gonboard = names(geneRanges)
    if (matchProbeNames) {
     pn = probeNames(mgr)
     gok = match(pn, gonboard, nomatch=0 )
     gok = gok[gok>0]
     geneRanges = geneRanges[ gok ]
     }
    if (length(gonboard)==0) stop("geneRanges must have non-null names")
    if (any(duplicated(gonboard))) stop("names of geneRanges must be unique")
#    if (is.null(ffind)) snm = unlist(lapply(mgr@fflist, rownames))
#    else snm = unlist(lapply(mgr@fflist[ffind], rownames))
    if (is(mgr, "eqtlTestsManager")) {
        snm = unlist(lapply(mgr@fflist, rownames))
     } else if (is(mgr, "cisTransDirector")) {
#       snm = unlist(snpIdList(mgrs(mgr)[[ffind]]))
       stop("this use case not covered with respect to ffind selection -- needs test")
     }

    iniRangeNames = names(snpRanges)
    onboard = intersect(snm, iniRangeNames)
    mm = match(onboard, iniRangeNames)
    snpRanges = snpRanges[mm]
    snps = containedSnps(geneRanges, snpRanges)
    gn = names(geneRanges)
    if (matchProbeNames) {
      ans = applier(1:length(snps), function(x) mgr[rsid(snps[[x]]), 
        probeId(gn[x])])
      names(ans) = names(snps)  # used to be names(geneRanges) ...
      } else {
      ans = applier(1:length(snps), function(x) mgr[rsid(snps[[x]]), ])
      names(ans) = names(snps)  # containedSnps propagated suitable range names
      }
    ans
}

