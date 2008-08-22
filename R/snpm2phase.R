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
}

parsePhPairs = function(fn) {
 li = readLines(fn)
 iinds = c(grep("IND", li), length(li)+1)
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

