
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
 

scoresInRanges = function (mgr, geneRanges, snpRanges, applier = lapply, ffind=NULL) 
{
 #
 # mgr is an eqtlTestsManager instance
 # geneRanges will typically be a GRanges extended for 'cis'
 # snpRanges can be a general snp metadata GRanges
 #  oct 22-- introduced ffind to reduce scope of mgr to be examined if geneRanges is limited to one chr
 #           also added filtering of geneRanges to genes available in mgr
 #
    pn = probeNames(mgr)
    gonboard = names(geneRanges)
    gok = match(pn, gonboard, nomatch=0 )
    gok = gok[gok>0]
    geneRanges = geneRanges[ gok ]
    if (is.null(ffind)) snm = unlist(lapply(mgr@fflist, rownames))
    else snm = unlist(lapply(mgr@fflist[ffind], rownames))
    iniRangeNames = names(snpRanges)
    onboard = intersect(snm, iniRangeNames)
    mm = match(onboard, iniRangeNames)
    snpRanges = snpRanges[mm]
    snps = containedSnps(geneRanges, snpRanges)
    gn = names(geneRanges)
    ans = applier(1:length(snps), function(x) mgr[rsid(snps[[x]]), 
        probeId(gn[x])])
    names(ans) = names(geneRanges)
    ans
}

