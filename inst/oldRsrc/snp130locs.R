
snp130locs = function(chr, start, end) {
 sess = browserSession()
 quer = ucscTableQuery(sess, "snp130", GenomicRanges(start, end, chr))
 tableName(quer) = "snp130"
 track(quer)
}

