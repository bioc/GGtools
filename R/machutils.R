

mrefhap2sm = function(gzfn, snpids) {
#
# use a reference haplotype data structure as supplied by MACH
# to create a SnpMatrix instance
#
 con = gzfile(gzfn)
 on.exit(close(con))
 s = scan(con, "")
 inds = seq(1,length(s),4) # assumes 2 records per id with one " " separating id and hap
 basicids = s[inds]
 c1 = s[inds+1]
 c2 = s[inds+3]
 c1vec = unlist(strsplit(c1, ""))
 c2vec = unlist(strsplit(c2, ""))
 idvec = rep(basicids,each=length(c1s[[1]]))
 snpidvec = rep(snpids, length(basicids))
 mat = cbind(idvec, snpidvec, c1vec, c2vec)
 targ = tempfile()
 on.exit(unlink(targ), add=TRUE)
 write.table(mat, file=targ, col.name=FALSE, row.names=FALSE, quote=FALSE)
 read.snps.long(targ, fields=c(sample=1, snp=2, allele1=3,
    allele2=4), codes="nucleotide", sample.id=basicids, snp.id=snpids)
 }
