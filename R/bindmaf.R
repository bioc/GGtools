
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

