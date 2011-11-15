geneRanges.annopk = function(ids, annopkg, extend=0) {
 require(annopkg, character.only=TRUE)
 anbase = gsub(".db", "", annopkg)
 chrmap = get(paste(anbase, "CHR", sep=""))
 stmap = get(paste(anbase, "CHRLOC", sep=""))
 endmap = get(paste(anbase, "CHRLOCEND", sep=""))
 chrs = mget(ids, chrmap)
 chrs = unlist(sapply(chrs, "[", 1))
 chrs = paste("chr", chrs, sep="")
 sts = mget(ids, stmap)
 sts = unlist(sapply(sts, "[", 1))
 ens = mget(ids, endmap)
 ens = unlist(sapply(ens, "[", 1))
 negsts = which(sts < 0)
 stra = rep("+", length(sts))
 stra[negsts] = "-"
 sts[negsts] = -sts[negsts]
 ens[negsts] = -ens[negsts]
 st = pmax(1,sts-extend)
 en = ens+extend
 st[is.na(sts)] = 1
 en[is.na(sts)] = 2
 #RangedData(IRanges(st,en), space=chrs, name=ids)
 tmp = GRanges(IRanges(st,en), seqnames=chrs, strand=factor(stra))
 names(tmp) = ids
 tmp
}
 
geneSyms = function(ids, annopkg) {
 require(annopkg, character.only=TRUE)
 anbase = gsub(".db", "", annopkg)
 symmap = get(paste(anbase, "SYMBOL", sep=""))
 unlist(sapply(mget(ids, symmap, ifnotfound=NA), "[", 1))
}

geneRanges = function(genomeOrPkgOrTxDb, chr=NULL, is.annopkg=FALSE, extend=0) {
 if (is.null(chr)) stop("must specify chr")
 if (is.annopkg) {
    if (!is(genomeOrPkgOrTxDb, "character")) stop("genomeOrPkgOrTxDb must be string if is.annopkg is TRUE")
    require(genomeOrPkgOrTxDb, character.only=TRUE)
    pkname = gsub(".db", "", genomeOrPkgOrTxDb)
    clenv = get(paste(pkname, "CHRLOC", sep=""))
    ids = mappedkeys(clenv)
    tmp = geneRanges.annopk(ids, genomeOrPkgOrTxDb, extend=extend)
    return(tmp[seqnames(tmp) == chr])
    }
 if (is(genomeOrPkgOrTxDb, "TranscriptDb")) txdb = genomeOrPkgOrTxDb
 else if (!(genomeOrPkgOrTxDb %in% c("hg18", "hg19")))
   stop("genomeOrPkgOrTxDb must be a TranscriptDb instance or %in% c('hg18', 'hg19')")
 else {
    txpk = paste("TxDb.Hsapiens.UCSC", genomeOrPkgOrTxDb, "knownGene",
      sep=".") 
    require(txpk, character.only=TRUE)
    txdb = get(txpk)
    }
 actseq = isActiveSeq(txdb)
 sn = names(actseq)
 if (!(chr %in% sn)) stop(paste("chr", chr, 
       " is not in names(isActiveSeq(txdb))"))
 actseq[] = FALSE
 actseq[chr] = TRUE
 isActiveSeq(txdb) = actseq
 txl = transcriptsBy(txdb, "gene")
 gn = names(txl)
 simpleExtents = function(r) IRanges(min(start(r)), max(end(r)))
 z = IRangesList(lapply(txl, function(z) simpleExtents(ranges(z))))
 strnd = sapply(txl, function(x)runValue(strand(x))[1])
 ans = IRanges::unlist(z)
 ans = GRanges(seqnames=chr, ans, strand=strnd)
 names(ans) = gn
 ans
}

.splitGR2GRL = function(g1, g2, ...) {
#
# a simple approach would be
#function(g1, g2) {
# nn = names(g2)
# ans = lapply(1:length(nn), function(z) subsetByOverlaps(g1, g2[z]))
# names(ans) = nn
# ans
#}
# but that is much slower
#
   if (is.null(names(g2))) names(g2) = as.character(1:length(g2))
   fo = findOverlaps(g1, g2, ...)
   mm = matchMatrix(fo)
   if (prod(dim(mm))==0) stop("no matchMatrix generated in findOverlaps")
#   fac = split(mm[,1], mm[,2])  # this reorders data so names of split list are ordered
   mmo = mm[order(mm[,2]), ]
   GN = names(g2)[mmo[,2]]
   fac = split(mmo[,1], GN) #mmo[,2])  # no further reordering
   fac = factor(rep(names(fac), sapply(fac,length)))
   g1 = g1[mmo[,1]]
   split(g1, fac)
}

 

