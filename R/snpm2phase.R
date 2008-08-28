setClass("phaseInput",
   representation(file4phase="character", exdir="character",
	sampid="character",
	snpid="character",
	locs="numeric"))

setMethod("show", "phaseInput", function(object) {
 cat("Input for PHASE:\n")
 cat(length(object@snpid), "SNP:\n")
 print(selectSome(object@snpid))
 cat(length(object@sampid), "Samples:\n")
 print(selectSome(object@sampid))
 cat("file for PHASE:\n")
 cat(object@file4phase, "\n")
 cat("[snpm2phase invoked in: ", object@exdir, "]\n")
})

setGeneric("invokePhase", function(x, cnum, parmstring, 
    globpname, where2run, doParse) standardGeneric("invokePhase"))

setMethod("invokePhase", c("phaseInput", "chrnum", "character", "character",
      "character", "logical"),
  function(x, cnum, parmstring, globpname, where2run, doParse) {
    if (!file.exists(globpname)) stop(paste(globpname, "does not exist but should be path for PHASE"))
    if (!file.exists(where2run)) stop(paste(where2run, "does not exist but should be folder where PHASE will run"))
    pinp = x@file4phase
    if (!file.exists(pinp)) stop(paste(pinp, "does not exist but should be input file for PHASE"))
    curwd = getwd()
    on.exit(setwd(curwd))
    setwd(where2run)
    path2 = function(x) gsub("[A-Za-z0-9]*$", "", x)
    execLine = paste(globpname, pinp, pout <- paste(pinp, "out", sep="."), parmstring)
    system(execLine)
    if (doParse) ans = parsePh.out(pout)
    ans
})

setMethod("invokePhase", c("snp.matrix", "chrnum", "character", "character",
      "character", "logical"),
  function(x, cnum, parmstring, globpname, where2run, doParse) {
    phin = snpm2phase(x, cnum, tempfile())
    if (!file.exists(globpname)) stop(paste(globpname, "does not exist but should be path for PHASE"))
    if (!file.exists(where2run)) stop(paste(where2run, "does not exist but should be folder where PHASE will run"))
    pinp = phin@file4phase
    if (!file.exists(pinp)) stop(paste(pinp, "does not exist but should be input file for PHASE"))
    curwd = getwd()
    on.exit(setwd(curwd))
    setwd(where2run)
    path2 = function(x) gsub("[A-Za-z0-9]*$", "", x)
    execLine = paste(globpname, pinp, pout <- paste(pinp, "out", sep="."), parmstring)
    system(execLine)
    if (doParse) ans = parsePh.out(pout)
    ans
})
    

snpm2phase = function(snpm, cnum, outfilename) {
 require(GGBase)
 con = file(outfilename,open="w")
 on.exit(close(con))
 snpid = colnames(snpm)
 sampid = rownames(snpm)
 nsnp = length(snpid)
 nsamp = length(sampid)
 loc = snpLocs.Hs(chrnum(cnum), rsid(snpid))["loc",]
 locl = paste("P", paste(loc, collapse=" "), sep = " ")
 sstr = paste(rep("S", nsnp), collapse="")
 writeLines(text=as.character(nsamp), con=con)
 writeLines(text=as.character(nsnp), con=con)
 writeLines(text=as.character(locl), con=con)
 writeLines(text=as.character(sstr), con=con)
 calls = as(snpm, "character")
 splc = list()
 for (i in 1:nsamp) {
    writeLines(text=as.character(sampid[i]), con=con)
    tmp = strsplit(calls[i,], "/")
    l1 = sapply(tmp, "[", 1)
    l2 = sapply(tmp, "[", 2)
    bad1 = which(nchar(l1) == 0 | is.na(l1))
    bad2 = which(nchar(l2) == 0 | is.na(l2))
    if (length(bad1)>0) l1[bad1] = "?"
    if (length(bad2)>0) l2[bad2] = "?"
    writeLines(text=paste(as.character(l1),collapse=""), con=con)
    writeLines(text=paste(as.character(l2),collapse=""), con=con)
 }
 cat(paste("wrote", outfilename), ".\n")
 new("phaseInput", snpid=snpid, sampid=sampid, locs=loc,
    file4phase=outfilename, exdir=getwd())
}

parsePhPairs = function(fn,subtok="IND") {
 li = readLines(fn)
 iinds = c(grep(subtok, li), length(li)+1)
 sep = rep(1:(length(iinds)-1), diff(iinds))
 strr = split(li, sep)
 ids = sapply(strr, "[", 1)
 data = lapply(strr, "[", -1) # drop id token
 names(data) = ids
 sdata = lapply(data, function(x) strsplit(x, " , "))
 probs = lapply(sdata, function(x) sapply(x, "[", 3))
 kp = sapply(probs, function(x)which.max(as.numeric(x)))
 tdata = list()
 for (i in 1:length(sdata))
   tdata[[i]] = sdata[[i]][[kp[i]]]
 names(tdata) = names(sdata)
 probs = as.numeric(sapply(tdata, "[", 3))
 names(probs) = names(sdata)
 tdata = lapply(tdata, "[", -3)
 list(tdata=tdata, probs=probs)
}

parsePh.out = function(fn) {
 li = readLines(fn)
 lsumini = grep("BEGIN LIST_SUMMARY", li)
 lsumend = grep("END LIST_SUMMARY", li)
 hlist = li[seq(lsumini+1, lsumend-1)]
 compressSp = function(x) gsub("[\ ]+"," ",x)
 dropLs = function(x) gsub("^ ", "", x)
 tmp = strsplit(dropLs(compressSp(hlist)), " ")
 inds = sapply(tmp, "[", 1)
 co = sapply(tmp, "[", 2)
 bpsumini = grep("BEGIN BESTPAIRS_SUMMARY", li)
 bpsumend = grep("END BESTPAIRS_SUMMARY", li)
 bpdat = li[seq(bpsumini+1, bpsumend-1)]
 id = sub(":.*", "", bpdat)
 tocom = sub(".*\\(", "", bpdat)
 tocom = sub(")", "", tocom)
 prs = strsplit(tocom, ",")
 el1 = as.numeric(sapply(prs, "[", 1))
 el2 = as.numeric(sapply(prs, "[", 2))
 ans = cbind( co[el1], co[el2])
 ans = t(apply(ans, 1, sort))
 al = list()
 for (i in 1:nrow(ans))
   al[[i]] = ans[i,]
 names(al) = id
 #ans = apply(ans, 1, paste, collapse=":")
 #names(ans) = id
 #ans
 pro = rep(NA, nrow(ans))
 names(pro) = id
 list(tdata=al, probs=pro)
}
