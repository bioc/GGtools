
 bindmaf = function(smpack="GGdata", smchr="20", obj) {
  rad = values(obj@scoregr)$radiusUsed[1]
  fr = fullreport(obj)
  pro = names(fr)
  values(fr)$probeid = pro
  togetv = values(fr)
  toget = togetv$snpid
  togetloc = togetv$snploc
  smls = getSS(smpack, smchr)
  probeanno = annotation(smls)
  glocenv = get(paste(gsub(".db", "", probeanno), "CHRLOC", sep=""))
  glocendenv = get(paste(gsub(".db", "", probeanno), "CHRLOCEND", sep=""))
  summ = col.summary(smList(smls)[[smchr]])
  rn = rownames(summ)
  ok = intersect(toget, rn)
  mafs = summ[ok,"MAF"]
  names(fr) = toget
  fr = fr[ok]  #  now we have the right set of probe ids
  okpr = values(fr)$probeid
  gstarts = sapply(mget(okpr, glocenv), "[", 1)
  gends = sapply(mget(okpr, glocendenv), "[", 1)
  strand = ifelse(gstarts<0, "-", "+")
  if (any(is.na(strand))) {
    warning("strand unknown for some gene, setting to +")
    strand[which(is.na(strand))] = "+" # FIXME
    }
  strand(fr) = strand
  values(fr)$MAF = mafs
#  will measure distance of SNP to midpoint of gene coding region
#  so that SNP within the gene have some distance too
  gmids = abs(gstarts) + .5*(abs(gends)-abs(gstarts))
  dist.mid = fr$snploc - gmids
  if (any(strand=="-"))   # use negative distance to denote SNP upstream (5') to midpoint of gene
    dist.mid[ which(strand == "-") ] = -dist.mid[which(strand=="-")]
  values(fr)$dist.mid = dist.mid
  fr
 }

richNull = function(..., MAFlb=.01, npc=10, radius=250000,
   nperm=1, innerFilt=function(x)x, outerFilt=function(x)x) {
  bigfilt = function(z) outerFilt(MAFfilter(clipPCs(permEx(innerFilt(z)), 1:npc), lower=MAFlb))
  inargs = list(...)
  if (any(names(inargs) %in% c("nperm", "npc", "radius", "MAFlb", "innerFilt"))) stop(
		"reserving argnames 'nperm', 'npc', 'radius', 'MAFlb', 'innerFilt', should not occur in ...")
  if (!(all(c("smpack", "chrnames") %in% names(inargs)))) stop("'smpack' and 'chrnames' are obligatory args")
  lapply(1:nperm, function(x)
    bindmaf(smpack=inargs$smpack,
            smchr=inargs$chrnames, 
            obj=best.cis.eQTLs( ..., smFilter=bigfilt, nperm=0, radius=radius )))
}
 

 meta.bindmaf = function(smpackvec=c("GGdata", "hmyriB36"), 
     smchr="20", obj) {
  rad = values(obj@scoregr)$radiusUsed[1]
  fr = fullreport(obj)
  pro = names(fr)
  values(fr)$probeid = pro
  togetv = values(fr)
  toget = togetv$snpid
  togetloc = togetv$snploc
  smls = getSS(smpackvec[1], smchr)
  probeanno = annotation(smls)
  glocenv = get(paste(gsub(".db", "", probeanno), "CHRLOC", sep=""))
  glocendenv = get(paste(gsub(".db", "", probeanno), "CHRLOCEND", sep=""))
#
# here need to generate minimum MAF over populations
#
  summ = unified.col.summary(smpackvec, smchr) #smList(smls)[[smchr]])
  rn = rownames(summ)
#
#
#
  ok = intersect(toget, rn)
  mafs = summ[ok,"MAF"]
  names(fr) = toget
  fr = fr[ok]  #  now we have the right set of probe ids
  okpr = values(fr)$probeid
  gstarts = sapply(mget(okpr, glocenv), "[", 1)
  gends = sapply(mget(okpr, glocendenv), "[", 1)
  strand = ifelse(gstarts<0, "-", "+")
  if (any(is.na(strand))) {
    warning("strand unknown for some gene, setting to +")
    strand[which(is.na(strand))] = "+" # FIXME
    }
  strand(fr) = strand
  values(fr)$MAF = mafs
#  will measure distance of SNP to midpoint of gene coding region
#  so that SNP within the gene have some distance too
  gmids = abs(gstarts) + .5*(abs(gends)-abs(gstarts))
  dist.mid = fr$snploc - gmids
  if (any(strand=="-"))   # use negative distance to denote SNP upstream (5') to midpoint of gene
    dist.mid[ which(strand == "-") ] = -dist.mid[which(strand=="-")]
  values(fr)$dist.mid = dist.mid
  fr
 }

unified.col.summary = function(smpackvec, smchr) {
 smats = lapply(smpackvec, function(x) getSS(x, smchr))
 npacks = length(smpackvec)
 oksn = colnames(smats[[1]])
 for (i in 2:npacks)
  oksn = intersect(oksn, colnames(smats[[i]]))
 summs = lapply(smats, function(x) col.summary(smList(x)[[1]]))
 outmafs = matrix(NA, nr=length(oksn), nc=npacks)
 for (ind in 1:npacks)  outmafs[,ind] = summs[[ind]][,"MAF"]
 ans = data.frame(MAF=apply(outmafs,1,min,na.rm=TRUE))
 rownames(ans) = oksn
 ans
}

meta.richNull = function(..., MAFlb=.01, npc=10, radius=250000,
   nperm=1, innerFilt=function(x)x, outerFilt=function(x)x) {
  bigfilt = function(z) outerFilt(MAFfilter(clipPCs(permEx(innerFilt(z)), 1:npc), lower=MAFlb))
  inargs = list(...)
  if (any(names(inargs) %in% c("nperm", "npc", "radius", "MAFlb", "innerFilt"))) stop(
		"reserving argnames 'nperm', 'npc', 'radius', 'MAFlb', 'innerFilt', should not occur in ...")
  if (!(all(c("smpack", "chrnames") %in% names(inargs)))) stop("'smpack' and 'chrnames' are obligatory args")
  lapply(1:nperm, function(x)
    bindmaf(smpack=inargs$smpack,
            smchr=inargs$chrnames, 
            obj=best.cis.eQTLs( ..., smFilter=bigfilt, nperm=0, radius=radius )))
}
