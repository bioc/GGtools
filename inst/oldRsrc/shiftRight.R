
snpcode = function(sms, chr, id)
  factor(as(snps(sms, chrnum(chr))[, id], "character"))

shiftRight = function( vbl, smlSet ) {
 vec = exprs(smlSet)[vbl,]
 pData(smlSet)[[vbl]] = vec
 smlSet
}

