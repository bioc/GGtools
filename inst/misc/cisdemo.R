
library(ceu1kg)
data(ceu1kgMeta_20)
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
gl = geneLimits( anno="lumiHumanAll.db" )

extendGR = function( gr, siz=1e6 )
  punion(shift(gr, -siz), shift(gr, siz), fill.gap=TRUE)

egl = extendGR(gl)

cisSnps = function(geneRanges, snpRanges) {
  mm = findOverlaps( geneRanges, snpRanges )@matchMatrix
  sinds = split( mm[,2], mm[,1] )  # split snp indexes by gene index
  ugi = unique(mm[,1])
  gn = names(geneRanges)[ugi]
  ans = lapply( 1:length(sinds), function(x) names(snpRanges)[sinds[[x]] ] )
  names(ans) = gn
  ans
}
