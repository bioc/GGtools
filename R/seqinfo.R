makeSeqinfo = function(genome="hg19") {
 if (!(genome=="hg19")) stop("only supporting hg19 at present")
 data(hg19.si.df)
 do.call(Seqinfo, hg19.si.df)
}
