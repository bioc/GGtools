
slimCTchr2SQLite = function(path,chr=1, dbname=paste("chr", chr, sep="")) {
 require(RSQLite, quietly=TRUE)
 if (!is.numeric(chr)) stop("chr must be numeric index (for Hs, in 1-24)")
 owd = getwd()
 on.exit(setwd(owd))
 tst = try(setwd(path))
 if (inherits(tst, "try-error")) stop("could not change to given folder")
 fn = dir(patt="\\+.rda$")
 obn = gsub(".rda", "", fn)
 load(fn[1])
 rsn = get(obn[1])$snpnames[[chr]]
 stats = get(obn[1])$chisq[[chr]]
 rm(list=obn[1])
 if (length(obn)>1) for (i in 2:length(obn))  {
   load(fn[i])
   stats = cbind(stats, get(obn[i])$chisq[[chr]])
   rm(list=obn[i])
   }
 stats = cbind(rsid=rsn,stats)
 dr = dbDriver("SQLite")
 con = dbConnect(dr, paste(dbname, ".db", sep=""))
 ans <- dbWriteTable(con, dbname, data.frame(stats))
 if (ans) cat("table written\n")
 NULL
}
 


slimCTchr2flat = function(path,chr=1, dbname=paste("chr", chr, sep="")) {
 if (!is.numeric(chr)) stop("chr must be numeric index (for Hs, in 1-24)")
 owd = getwd()
 on.exit(setwd(owd))
 tst = try(setwd(path))
 if (inherits(tst, "try-error")) stop("could not change to given folder")
 fn = dir(patt="\\+.rda$")
 obn = gsub(".rda", "", fn)
 load(fn[1])
 rsn = get(obn[1])$snpnames[[chr]]
 stats = get(obn[1])$chisq[[chr]]
 rm(list=obn[1])
 if (length(obn)>1) for (i in 2:length(obn))  {
   load(fn[i])
   stats = cbind(stats, get(obn[i])$chisq[[chr]])
   rm(list=obn[i])
   }
 stats = cbind(rsid=rsn,stats)
 ofile = gzfile(paste(dbname, ".txt.gz", sep=""))
 write.table(stats, ofile)
 cat("table written\n")
 NULL
}
 

