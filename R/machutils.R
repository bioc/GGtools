

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
 idvec = rep(basicids,each=length(snpids))
 snpidvec = rep(snpids, length(basicids))
 mat = cbind(idvec, snpidvec, c1vec, c2vec)
 targ = tempfile()
 on.exit(unlink(targ), add=TRUE)
 write.table(mat, file=targ, col.name=FALSE, row.names=FALSE, quote=FALSE)
 read.snps.long(targ, fields=c(sample=1, snp=2, allele1=3,
    allele2=4), codes="nucleotide", sample.id=basicids, snp.id=snpids)
 }


sm2ped = function(sm, snpsupp, missing.code="N", 
  family, person, father, mother, sex ) { 
  tmp = whmu( sm, snpsupp, missing.code=missing.code, transpose=TRUE )
  if (missing(family)) family=1:nrow(sm)
  if (missing(person)) person=1:nrow(sm)
  if (missing(father)) father=rep(0,nrow(sm))
  if (missing(mother)) mother=rep(0,nrow(sm))
  if (missing(sex)) sex=rep(1,nrow(sm)) # female=1, male=2
  cbind(family, person, father, mother, sex, tmp)
}
  
# transform unphased hapmap data from snpMatrix + support to char matrix with one row per subject
whmu = function( sm, sup, missing.code="-", markerdf=NULL, transpose=FALSE) {
# to generate beagle unphased input write.table(OO, file="OO.txt", quote=FALSE, row.names=FALSE)
 snames = colnames(sm)
 ids = rownames(sm)
 nsnp = ncol(sm)
 nsup = nrow(sup)
 if (nsnp != nsup) stop("number of snps in sm not equal to number of support lines")
 nsub = nrow(sm)
 ass = as.character(sup[,"Assignment"])
 calls = strsplit(ass, "/")
 calls = sapply(calls, function(x) {x[!x%in%c("A", "C", "G", "T")] = missing.code; x})
 smn = as(sm, "numeric")  # 0,1,2 and NA
 a2 = a1 = smn
 a1[a1 <=  1] = 1 # 0 means homozygous early code, 1 means het
 a2[smn == 0] = 1
 a2[smn == 1] = 2
 r1 = sapply(1:nsub, function(x) {
    calls[ cbind(a1[x,], 1:nsnp) ]
    })
 r2 = sapply(1:nsub, function(x) {
    calls[ cbind(a2[x,], 1:nsnp) ]
    })
 colnames(calls) = snames
 if (any(is.na(r1))) r1[is.na(r1)] = missing.code
 if (any(is.na(r2))) r2[is.na(r2)] = missing.code
 rownames(r1) = snames
 colnames(r1) = ids
 rownames(r2) = snames
 colnames(r2) = ids
 if (transpose) {
  r1 = t(r1)
  r2 = t(r2)
 }
 final = matrix(rbind(r1,r2), nr=nrow(r1))
 if (transpose) {
   colnames(final) = rep(snames,each=2)
   return(final)
   }
 colnames(final) = rep(ids, each=2)
 rownames(final) = snames
 ncalls = calls
if (!is.null(markerdf)) {
 okid = intersect(snames, as.character(markerdf[,1]))
 final = final[okid,]
 ncalls = calls[,okid]
 rownames(markerdf) = as.character(markerdf[,1])
 mdf = markerdf[okid,]
 mmm = apply(mdf[,3:4], 1, paste, collapse="")
 nnn = apply(ncalls,2,paste,collapse="")
 ok = which(mmm == nnn)
 if (length(ok) == 0) stop("no matches to marker file")
 final = final[ok,]
 }
 data.frame(I="M", id=rownames(final), final, check.names=FALSE)
}
 
