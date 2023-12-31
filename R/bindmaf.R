
 bindmaf.legacy = function(smpack="GGdata", smchr="20", obj, SSgen=GGBase::getSS) {
#
  rad = values(obj@scoregr)$radiusUsed[1]
  fr = fullreport(obj)
  fr = fr[ which(as.character(seqnames(fr)) == smchr) ]
  pro = names(fr)
  values(fr)$probeid = pro
  togetv = values(fr)
  toget = togetv$snpid
  togetloc = togetv$snploc
  smls = SSgen(smpack, smchr)
  probeanno = annotation(smls)
  require(probeanno, character.only=TRUE)
  glocenv = get(paste(gsub(".db", "", probeanno), "CHRLOC", sep=""))
  glocendenv = get(paste(gsub(".db", "", probeanno), "CHRLOCEND", sep=""))
  summ = col.summary(smList(smls)[[smchr]])
  rn = rownames(summ)
#  ok = toget[ toget %in% rn ] # intersect(toget, rn) -- retain shared SNP
  if (!all(toget %in% rn)) stop("some SNP not available in SSgen result ... shouldn't happen")
#  names(fr) = toget
#  fr = fr[ok]  #  now we have the right set of probe ids
  if (!all.equal(names(fr), values(fr)$probeid)) stop("probeides went out of sync")
  mafs = summ[match(toget, rn),"MAF"]
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
  mindist = ifelse(fr$snploc >= abs(gstarts) &
      fr$snploc <= abs(gends), 0, pmin(fr$snploc -abs(gstarts), fr$snploc-abs(gends)))
  dist.mid = fr$snploc - gmids
  if (any(strand=="-"))   # use negative distance to denote SNP upstream (5') to midpoint of gene
    dist.mid[ which(strand == "-") ] = -dist.mid[which(strand=="-")]
  values(fr)$dist.mid = dist.mid
  values(fr)$mindist = mindist
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
     smchr="20", obj, usemaxMAF=FALSE, SSgen=GGBase::getSS) {
  rad = values(obj@scoregr)$radiusUsed[1]
  fr = fullreport(obj)
  fr = fr[ which(as.character(seqnames(fr)) == smchr) ]
  pro = names(fr)
  values(fr)$probeid = pro
  togetv = values(fr)
  toget = togetv$snpid  # these need not be unique.  two genes can 
           # share one best snp
  togetloc = togetv$snploc
  smls = SSgen(smpackvec[1], smchr)
  probeanno = annotation(smls)
  require(probeanno, character.only=TRUE)
  glocenv = get(paste(gsub(".db", "", probeanno), "CHRLOC", sep=""))
  glocendenv = get(paste(gsub(".db", "", probeanno), "CHRLOCEND", sep=""))
#
# here need to generate minimum MAF over populations
#
  summ = unified.col.summary(smpackvec, smchr, usemax=usemaxMAF) #smList(smls)[[smchr]])
  rn = rownames(summ)
#
#
#
#  ok = toget[ toget %in% rn ] # intersect(toget, rn) -- retain shared SNP
  if (!all(toget %in% rn)) stop("some SNP not available in SSgen result ... shouldn't happen")
  mafs = summ[match(toget, rn),"MAF"]
#  names(fr) = toget
#  fr = fr[ok]  #  now we have the right set of probe ids
  if (!all.equal(names(fr), values(fr)$probeid)) stop("probeides went out of sync")
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

unified.col.summary = function(smpackvec, smchr, usemax=FALSE,
   SSgen=GGBase::getSS) {
##
## usemax = TRUE -> take max MAF over all populations
## otherwise compute the MAF for the combined samples
##
 smats = lapply(smpackvec, function(x) smList(SSgen(x, smchr))[[1]])
 npacks = length(smpackvec)
 oksn = colnames(smats[[1]])
 for (i in 2:npacks)
  oksn = intersect(oksn, colnames(smats[[i]]))
 for (i in 1:npacks)
  smats[[i]]  = smats[[i]][,oksn]
 if (usemax) {
   summs = lapply(smats, col.summary)
   outmafs = matrix(NA, nrow=length(oksn), ncol=npacks)
   for (ind in 1:npacks)  outmafs[,ind] = summs[[ind]][oksn,"MAF"]
   ans = data.frame(MAF=apply(outmafs,1,max,na.rm=TRUE))
   }
 else {
   fullmat = do.call( BiocGenerics::rbind, smats )
   summs = col.summary( fullmat )
   ans = data.frame(MAF=summs[oksn, "MAF"])
   }
 rownames(ans) = oksn
 ans
}

meta.richNull = function(..., MAFlb=.01, npc=10, radius=250000,
   nperm=1, innerFilt=function(x)x, outerFilt=function(x)x) {
  bigfilt = function(z) outerFilt(MAFfilter(clipPCs(permEx(innerFilt(z)), 1:npc), lower=MAFlb))
  inargs = list(...)
  if (any(names(inargs) %in% c("nperm", "npc", "radius", "MAFlb", "innerFilt"))) stop(
		"reserving argnames 'nperm', 'npc', 'radius', 'MAFlb', 'innerFilt', should not occur in ...")
  if (!(all(c("smpackvec", "chrnames") %in% names(inargs)))) stop("'smpack' and 'chrnames' are obligatory args")
  bigfiltList = lapply(1:length(inargs$smpackvec), function(x) bigfilt)
  lapply(1:nperm, function(x)
    meta.bindmaf(smpackvec=inargs$smpackvec,
            smchr=inargs$chrnames, 
            obj=meta.best.cis.eQTLs( ..., SMFilterList=bigfiltList, nperm=0, radius=radius )))
}


bindmaf.simple = function(smpack, smchr, fr, SSgen=GGBase::getSS, rad, conf) {
# at this point, the object will have Homo.sapiens hg19 seqinfo
# but smchr may use a different vocabulary
#
# we also use the extraProps component of the configuration
# object to add additional metadata
#
   smchr.init = smchr  # ???
  smchr = gsub("chr", "", smchr)
  smchr = paste0("chr", smchr)
    fr = fr[which(as.character(seqnames(fr)) == smchr)]
    pro = names(fr)
    values(fr)$probeid = pro
    togetv = values(fr)
    toget = togetv$snp
    togetloc = togetv$snplocs
    smtoget = paste0(smchrpref(conf), chrnames(conf)[1])
    smls = SSgen(smpack, smchr.init, smFilter=smFilter(conf),
      exFilter=exFilter(conf))
    probeanno = annotation(smls)
    require(probeanno, character.only = TRUE)
    glocenv = get(paste(gsub(".db", "", probeanno), "CHRLOC", 
        sep = ""))
    glocendenv = get(paste(gsub(".db", "", probeanno), "CHRLOCEND", 
        sep = ""))
    summ = col.summary(smList(smls)[[1]])
    rn = rownames(summ)
    if (!all(toget %in% rn)) 
        stop("some SNP not available in SSgen result ... shouldn't happen")
    if (!all.equal(names(fr), values(fr)$probeid)) 
        stop("probeides went out of sync")
    mafs = summ[match(toget, rn), "MAF"]
    okpr = values(fr)$probeid
    gstarts = sapply(mget(okpr, glocenv), "[", 1)
    gends = sapply(mget(okpr, glocendenv), "[", 1)
    strand = ifelse(gstarts < 0, "-", "+")
    if (any(is.na(strand))) {
        warning("strand unknown for some gene, setting to +")
        strand[which(is.na(strand))] = "+"
    }
    strand(fr) = strand
    values(fr)$MAF = mafs
    gmids = abs(gstarts) + 0.5 * (abs(gends) - abs(gstarts))
    dist.mid = fr$snplocs - gmids
    mindist = ifelse(fr$snplocs >= abs(gstarts) & 
               fr$snplocs <= abs(gends), 
                   0, 
                   pmin(abs(fr$snplocs-abs(gstarts)), abs(fr$snplocs-abs(gends))))
    if (any(strand == "-")) 
        dist.mid[which(strand == "-")] = -dist.mid[which(strand == 
            "-")]
    values(fr)$dist.mid = dist.mid
    values(fr)$mindist = mindist
    values(fr)$genestart = gstarts
    values(fr)$geneend = gends
    fr
}


bindprops = function( config, fr ) {
#   smpack, smchr, fr, SSgen=GGBase::getSS, rad, conf) {
# at this point, the object will have Homo.sapiens hg19 seqinfo
# but smchr may use a different vocabulary
#
# we also use the extraProps component of the configuration
# object to add additional metadata
#
  smpack = smpack(config)
  smchr = chrnames(config)
  SSgen = SSgen(config)
  rad = radius(config)
   smchr.init = smchr
  smchr = gsub("chr", "", smchr)
  smchr = paste0("chr", smchr)
    fr = fr[which(as.character(seqnames(fr)) == smchr)]
    pro = names(fr)
    values(fr)$probeid = pro
    togetv = values(fr)
    toget = togetv$snp
    togetloc = togetv$snplocs
    togetsm = paste0(smchrpref(config), chrnames(config)[1])
    smls = SSgen(smpack, togetsm, exFilter=exFilter(config))
    probeanno = annotation(smls)
    require(probeanno, character.only = TRUE)
    glocenv = get(paste(gsub(".db", "", probeanno), "CHRLOC", 
        sep = ""))
    glocendenv = get(paste(gsub(".db", "", probeanno), "CHRLOCEND", 
        sep = ""))
    summ = col.summary(smList(smls)[[1]])
    rn = rownames(summ)
    if (!all(toget %in% rn)) 
        stop("some SNP not available in SSgen result ... shouldn't happen")
    if (!all.equal(names(fr), values(fr)$probeid)) 
        stop("probeides went out of sync")
    mafs = summ[match(toget, rn), "MAF"]
    okpr = values(fr)$probeid
    gstarts = sapply(mget(okpr, glocenv), "[", 1)
    gends = sapply(mget(okpr, glocendenv), "[", 1)
    strand = ifelse(gstarts < 0, "-", "+")
    if (any(is.na(strand))) {
        warning("strand unknown for some gene, setting to +")
        strand[which(is.na(strand))] = "+"
    }
    strand(fr) = strand
    values(fr)$MAF = mafs
    gmids = abs(gstarts) + 0.5 * (abs(gends) - abs(gstarts))
    dist.mid = fr$snplocs - gmids
    mindist = ifelse(fr$snplocs >= abs(gstarts) & 
               fr$snplocs <= abs(gends), 
                   0, 
                   pmin(abs(fr$snplocs-abs(gstarts)), abs(fr$snplocs-abs(gends))))
    if (any(strand == "-")) 
        dist.mid[which(strand == "-")] = -dist.mid[which(strand == 
            "-")]
    values(fr)$dist.mid = dist.mid
    values(fr)$mindist = mindist
    values(fr)$genestart = gstarts
    values(fr)$geneend = gends
    extraProps(config)(fr)
}
