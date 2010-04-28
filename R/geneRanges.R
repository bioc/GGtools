geneRanges = function(ids, annopkg, extend=0) {
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
 sts[negsts] = -sts[negsts]
 ens[negsts] = -ens[negsts]
 st = pmax(1,sts-extend)
 en = ens+extend
 st[is.na(sts)] = 1
 en[is.na(sts)] = 2
 RangedData(IRanges(st,en), space=chrs, name=ids)
}
 
geneSyms = function(ids, annopkg) {
 require(annopkg, character.only=TRUE)
 anbase = gsub(".db", "", annopkg)
 symmap = get(paste(anbase, "SYMBOL", sep=""))
 unlist(sapply(mget(ids, symmap, ifnotfound=NA), "[", 1))
}
 
