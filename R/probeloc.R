
# futuristic use of self-reference for a feature, for performance analysis

probeGen = setRefClass("probe", 
   fields=list(probename="character", annodbname="character"), 
   methods=list(
     location=function() {
		require(annodbname, character.only=TRUE)
		stem = gsub(".db", "", annodbname)
		stenv = get(paste(stem, "CHRLOC", sep=""))
		enenv = get(paste(stem, "CHRLOCEND", sep=""))
		st= get(probename, stenv)
		en= get(probename, enenv)
		chr = names(st)[1]
		GRanges(seqnames=chr, IRanges(st,en))
		},
     chr = function() {
		require(annodbname, character.only=TRUE)
		stem = gsub(".db", "", annodbname)
		stenv = get(paste(stem, "CHR", sep=""))
		get(probename, stenv)[1]
		}
))

# bread and butter location support
		
probeChromosomes = function(sms) {
	fn = featureNames(sms)
	stem = gsub(".db", "", annotation(sms))
	stenv = get(paste(stem, "CHR", sep=""))
	chrs = mget(fn, stenv, ifnotfound=NA)
	sapply(chrs, "[", 1)
}

probeLocations = function(sms, extend=0) {
	fn = featureNames(sms)
	stem = gsub(".db", "", annotation(sms))
	stenv = get(paste(stem, "CHRLOC", sep=""))
	st = mget(fn, stenv, ifnotfound=NA)
	st1 = sapply(st, "[", 1)
	bad = which(is.na(st1))
	ch1 =  lapply(st, function(x)names(x)[1])  # will be list for later subsetting
	enenv = get(paste(stem, "CHRLOCEND", sep=""))
	en = mget(fn, enenv, ifnotfound=NA)
	en1 = sapply(en, "[", 1)
	bad = union(bad, which(is.na(en1)))
	if (any(bad)) {
		warning(paste("dropped ", length(bad), " with no location"))
		st1 = st1[-bad]
		en1 = en1[-bad]
		ch1 = unlist(ch1[-bad])
	} else ch1 = unlist(ch1)
	strand = ifelse(st1<0, "-", "+")
        vec2rle = function(x) Rle(factor(x), rep(1, length(x)))
	ans = GRanges(seqnames=vec2rle(ch1), IRanges(abs(st1), abs(en1))+extend, strand=vec2rle(strand))
	if (length(bad) > 0) names(ans) = fn[-bad]
   	else names(ans) = fn
	ans
}

probeSequences = function(sms) {
	fn = featureNames(sms)
	stem = gsub(".db", "", annotation(sms))
	stenv = try(get(lkps <- paste(stem, "PROBESEQUENCE", sep="")))
	if (inherits(stenv, "try-error")) stop(paste(lkps, "not available ... only for illuminaHumanv...?"))
	mget(fn, stenv, ifnotfound=NA)
}

snpLocations = function(sms, snpLocGRanges, grsnpid = "RefSNP_id" ) {
	slc = names(smList(sms))
	sqn = seqlevels(snpLocGRanges)
	if (!all(slc %in% sqn)) stop(paste(
              "names(smList(sms)) [",
                  paste(sQuote(names(smList(sms))), collapse=","), 
                  "]  returns values not found in seqlevels(snpLocGRanges).  Please reconcile."))
	locl = list()
	for (i in 1:length(slc)) {
		tsn = colnames(smList(sms)[[i]])
		tsn = gsub("rs", "", tsn)
		ma = na.omit(match(tsn, elementMetadata(snpLocGRanges)[[grsnpid]]))
		locl[[i]] = snpLocGRanges[ma,]
		}
	names(locl) = slc
	locl
}

proximityList = function(sms, smlind=1, snpLocGRanges, grsnpid = "RefSNP_id", probeLocExtend=0,
   glocTransform = function(x)x) {
   gloc = glocTransform(probeLocations(sms,extend=probeLocExtend))
   sloc = snpLocations(sms=sms, snpLocGRanges=snpLocGRanges, grsnpid=grsnpid)[[smlind]]
   snames = paste("rs", elementMetadata(sloc)[[grsnpid]], sep="")
   ov = matchMatrix(findOverlaps(gloc, sloc))
   snpsmatched = snames[ov[,2]]
   probesmatched = names(gloc)[ov[,1]]
   split(snpsmatched,probesmatched)
}


getGene2SnpList = function(sms, chr, genome, radius=50000,
   additionalSNPGR=NULL, useTxDb=FALSE) {
#
# will return list named with probe names; each probe has a vector
#   of associated SNP names
# based on the bioconductor annotations
# note that "imputed" SNP in the sms may need additionalSNPGR info
#
# comments: return GRanges for genes and then a GRangesList for snp
# consider a generalization of GRangesList so that queries can be
# against list node elements
#
# first acquire the SNP locations
#
 if (length(chr)>1) stop("chr must be scalar")
 chr = as.character(chr)
 if (!(chr %in% as.character(1:22))) stop("chr must be in as.character(1:22)")
 if (genome == "hg18") spack = "SNPlocs.Hsapiens.dbSNP.20090506"
 else if (genome == "hg19") spack = "SNPlocs.Hsapiens.dbSNP.20100427"
 else stop("only using hg18 and hg19 to select SNPlocs")
 require(spack, character.only=TRUE)
 if (packageVersion(spack) < package_version("0.99.6")) stop(
	"please install SNPlocs package with version at least 0.99.6 [these have GRanges conversion built in]")
 seqn = do.call(":::", list(spack, "SEQNAMES"))  # 
# here we will assume that prefix of seqn is either "ch" or "chr"
 nseqn = gsub("ch", "", gsub("chr", "", seqn))
 seqnToUseInd = which(nseqn == chr)
 if (length(seqnToUseInd)==0) stop(paste("no match of chr submitted in call to SEQNAMES of", spack))
 seqnToUse = seqn[seqnToUseInd]
 getter = do.call("::", list(spack, "getSNPlocs"))  # for multiple SNPlocs on searchlist
 snpgr = getter( seqnToUse, as.GRanges=TRUE )  # LAST USE of seqnToUse, reset to chr prefix below
 elementMetadata(snpgr)$RefSNP_id = paste("rs", elementMetadata(snpgr)$RefSNP_id, sep="")
#
# deal with possible additional snpGR for imputed loci
#
# force chr prefix
 nnsl = gsub("[^0-9]", "", seqlevels(snpgr))  # numeric only
 seqlevels(snpgr) = make.names(paste("chr", nnsl, sep=""),unique=TRUE)
 if (!is.null(additionalSNPGR)) {
   emn = names(elementMetadata(additionalSNPGR))
   if (!all.equal(emn, names(elementMetadata(snpgr)))) stop(paste("need elementMetadata on additionalSNPGR to have columns", paste(names(elementMetadata(snpgr)), collapse=", ")))
   snpgr = c(snpgr, additionalSNPGR)
 }
#
#  deal with gene location
#
#  first determine chromosome residence
#
 ganno = annotation(sms)
 require(ganno, character.only=TRUE)
 mapper = get(paste(gsub(".db", "", ganno), "CHR", sep=""))
 ponc = get(chr, revmap(mapper))
 ponc = intersect(featureNames(sms), ponc)
 if (length(ponc) == 0) stop(paste("sms contains no probes on chromosome ", chr))
#
# now obtain coordinates
#
 if (!useTxDb) {  # will use bioc chip annotation location data -- typically hg19 for
                  # current R  
   if (!(genome == "hg19")) warning("you are using bioconductor locations; current are hg19-based")
   locmapper = get(paste(gsub(".db", "", ganno), "CHRLOC", sep=""))
   lenmapper = get(paste(gsub(".db", "", ganno), "CHRLENGTHS", sep=""))
   locendmapper = get(paste(gsub(".db", "", ganno), "CHRLOCEND", sep=""))
   pstart = abs(sapply(mget(ponc, locmapper, ifnotfound=NA), "[", 1))
   pend = abs(sapply(mget(ponc, locendmapper, ifnotfound=NA), "[", 1))
   if (any(is.na(c(pstart,pend)))) {
     bad = union( which(is.na(pstart)), which(is.na(pend)))
     pstart = pstart[-bad]
     pend = pend[-bad]
     ponc = ponc[-bad]
     }
   genelocs = GRanges(seqnames=paste("chr", chr, sep=""), IRanges(pstart,pend))
   maxlen = lenmapper[chr]
   } else { 
#
# using TxDb 
#
   egmapper = get(paste(gsub(".db", "", ganno), "ENTREZID", sep=""))
   egids = sapply(mget(ponc, egmapper, ifnotfound=NA), "[", 1)
   if (any(is.na(egids))) {
      badp = which(is.na(egids))
      ponc = ponc[-badp]
      egids = egids[-badp]
      }
   dupeg = which(duplicated(egids))
   if (length(dupeg)>0) egids = egids[-dupeg]
   TxDb2get = paste("TxDb.Hsapiens.UCSC", genome, "knownGene", sep=".")
   require(TxDb2get, character.only=TRUE)
   dbobj = get(TxDb2get)
   #txseqn = names(isActiveSeq(dbobj))
   #toset = (txseqn == paste("chr", chr, sep=""))  # crucial filtering
   #names(toset) = txseqn
   #isActiveSeq(dbobj) = toset
setSoloSeq = function(tx, ind) {
   curactive = isActiveSeq(tx)
   kpns = names(curactive)
   curactive[] = FALSE
   curactive[ind] = TRUE
   names(curactive) = kpns
   isActiveSeq(tx) = curactive
   tx
}
   dbobj = setSoloSeq(dbobj, paste("chr", chr, sep=""))
   tx = transcriptsBy(dbobj, "gene")
   egwaddr = names(tx)
   badeg = which(!(egids %in% egwaddr))
   if (length(badeg)>0) {
      egids = egids[-badeg]
      ponc = ponc[-badeg]
      }
   txok = reduce(tx[egids])
   alll = elementLengths(txok)
   if (any(alll != 1)) {
     warning("some entrez genes selected in TxDb had disjoint transcripts; picking first for addressing")
     txok = txok[-which(alll != 1)]
     }
   egnames = names(txok)
   genelocs = unlist(txok) #unlist(GRangesList(txok))
#   expand to probes
   pids = try(mget(egnames, revmap(egmapper) ))  # should not have NA
   if (inherits(pids, "try-error")) stop("unexpected NA in EG lookup, bug.")
   repv = sapply(pids, length)
   ponc = unlist(pids)
   genelocs = rep(genelocs,repv)
   strand(genelocs) = "*"
   maxlen = seqlengths(txok)
   }
#
# add radius extension
 exlocs = ranges(genelocs) + radius
 hasnonp = which(start(exlocs)<1)
 if (length(hasnonp) > 0) start(exlocs[hasnonp]) = 1
 hasoverhang = which(end(exlocs)>maxlen)
 if (length(hasoverhang) > 0) send(exlocs[hasoverhang]) = maxlen
 ranges(genelocs) = exlocs
 ol = findOverlaps(genelocs, snpgr)
 mmol = matchMatrix(ol)
 g2use = unique(mmol[,1])
 ponc = ponc[g2use]
 snpinds = split(mmol[,2], mmol[,1])
 snid = elementMetadata(snpgr)$RefSNP_id
 ans = lapply(snpinds, function(x) snid[x])
 names(ans) = ponc
 ans
}

 
