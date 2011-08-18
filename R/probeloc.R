
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

probeLocations = function(sms) {
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
	}
	strand = ifelse(st1<0, "-", "+")
	ans = GRanges(seqnames=ch1, IRanges(abs(st1), abs(en1)), strand=strand)
	names(ans) = fn[-bad]
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
	if (!all(slc %in% sqn)) stop("names(smList(sms)) returns values not found in seqlevels(snpLocGRanges).  Please reconcile.")
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
