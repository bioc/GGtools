filterVCF = function(gzpath, chrom, nrec=NULL, outfile=NULL, return.pipe=TRUE, list.chromnames=FALSE, tabixcmd = "tabix") {
#
# purpose is to simplify access to records in VCF via tabix index
#   can return a pipe with parameterized call to tabix
#   can direct tabix output to a file
#   can get tabix output internally as character vector
#
 if (list.chromnames) return(system(paste(tabixcmd, "-l",  gzpath), intern=TRUE))
 if (missing(chrom)) stop("must specify a chromosome")
 if (return.pipe) return(pipe(paste(tabixcmd, gzpath, chrom), "r"))
 redirec = ""
 recfilt = ""
 if (!is.null(nrec)) recfilt = paste("|head -", nrec, sep="")
 if (!is.null(outfile)) redirec = paste(">", outfile)
 system(paste(tabixcmd, gzpath, chrom, recfilt, redirec), intern=TRUE)
}
setClass("metaVCF", contains="character")

getMetaVCF = function (gzcon, maxnlines = 20, final = "^#CHROM") 
{
    if (!is(gzcon, "gzfile")) 
        stop("gzcon must inherit from gzfile; see ?gzfile")
    if (summary(gzcon)$opened != "opened") 
        stop("please pass an opened connection")
    meta = NULL
    for (i in 1:maxnlines) {
        meta = c(meta, tmp <- readLines(gzcon, n = 1))
        if (length(grep(final, tmp)) > 0) 
            break
    }
    new("metaVCF", meta)
}

setMethod("show", "metaVCF", function(object) cat(object, sep="\n"))


sampleIDs = function(metavec, ndrop=9) 
    strsplit( metavec[length(metavec)], "\t" )[[1]][-c(1:ndrop)]

#setMethod("sampleNames", c("metaVCF", "missing"), function(object, ...)
#  .sampleNames( object, ndrop=9 ) )
parseVCFrec = function(rec, nmetacol=9 ) {
 vec = strsplit(rec, "\t")[[1]]
 meta = vec[1:nmetacol]
 calls = vec[-c(1:nmetacol)]
 nalt = strsplit(calls, "")
 nums = lapply(nalt, "[", c(1,3))
 nalt = sapply(nums, function(x) 2-sum(x=="0"))
 if (any(is.na(nalt))) nalt[which(is.na(nalt))]=-1
 nalt = nalt+1
 chr = meta[1]
 id = meta[3]
 loc = meta[2]
 if (id == ".") id = paste("chr", chr, ":", loc, sep="")
 list(chr=chr, id=id, loc=loc, ref=meta[4], alt=meta[5], depth=meta[8],
   calls=as.raw(nalt))
}

vcf2sm = function(gzpath, chrom, tabixcmd = "tabix", nmetacol=9, verbose=FALSE) {
 require(snpMatrix)
 mm = getMetaVCF( gzfile(gzpath, "r") )
 sampids = sampleIDs(mm, ndrop=nmetacol)
 on.exit(close(fpipe))
 fpipe = filterVCF( gzpath, chrom, return.pipe = TRUE, tabixcmd = tabixcmd )
 out = list()
 i = 1
 while ( length(tmp <- readLines(fpipe, n=1))>0) {
  out[[i]] = parseVCFrec( tmp, nmetacol=nmetacol )
  if (verbose & (i%%1000)==0) cat(i)
  i = i+1
  }
 rsid = sapply(out, "[[", "id")
 nsnp = length(out)
 mat = matrix(as.raw(0), nr=length(sampids), ncol=nsnp)
 for (i in 1:nsnp) mat[,i] = out[[i]]$calls
 rownames(mat) = sampids
 colnames(mat) = rsid
 new("snp.matrix", mat)
}