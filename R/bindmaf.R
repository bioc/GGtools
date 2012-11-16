
 bindmaf = function(smpack="GGdata", smchr="20", obj) {
  rad = values(obj@scoregr)$radiusUsed[1]
  fr = fullreport(obj)
  pro = names(fr)
  values(fr)$probeid = pro
  togetv = values(fr)
  toget = togetv$snpid
  togetloc = togetv$snploc
  smls = getSS(smpack, smchr)
  summ = col.summary(smList(smls)[[smchr]])
  rn = rownames(summ)
  ok = intersect(toget, rn)
  mafs = summ[ok,"MAF"]
  names(fr) = toget
  fr = fr[ok]
  values(fr)$MAF = mafs
  ranges(fr) = ranges(fr)-rad
  mindist = pmin(abs(fr$snploc-start(fr)), abs(fr$snploc-end(fr)))
  swithing = which((fr$snploc >= start(fr)) & (fr$snploc <= end(fr)))
  if (length(swithing)>0) mindist[swithing] = 0
  values(fr)$dist.out = mindist
  fr
 }

richNull = function(..., MAFlb=.01, npc=10, radius=250000,
   nperm=1, innerFilt=function(x)x) {
  bigfilt = function(z) MAFfilter(clipPCs(permEx(innerFilt(z)), 1:npc), lower=MAFlb)
  inargs = list(...)
  if (any(names(inargs) %in% c("nperm", "npc", "radius", "MAFlb", "innerFilt"))) stop(
		"reserving argnames 'nperm', 'npc', 'radius', 'MAFlb', 'innerFilt', please resubmit without using these")
  if (!(all(c("smpack", "chrnames") %in% names(inargs)))) stop("'smpack' and 'chrnames' are obligatory args")
  lapply(1:nperm, function(x)
    bindmaf(smpack=inargs$smpack,
            smchr=inargs$chrnames, 
            obj=best.cis.eQTLs( ..., smFilter=bigfilt, nperm=1 )))
}
 
