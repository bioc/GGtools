#filterVCF = function(gzpath, chrom, nrec=NULL, outfile=NULL, return.pipe=TRUE, list.chromnames=FALSE, tabixcmd = "tabix") {
##
## purpose is to simplify access to records in VCF via tabix index
##   can return a pipe with parameterized call to tabix
##   can direct tabix output to a file
##   can get tabix output internally as character vector
##
# if (list.chromnames) return(seqnamesTabix(gzpath))
# if (missing(chrom)) stop("must specify a chromosome")
# if (return.pipe) return(pipe(paste(tabixcmd, gzpath, chrom), "r"))
# redirec = ""
# recfilt = ""
# if (!is.null(nrec)) recfilt = paste("|head -", nrec, sep="")
# if (!is.null(outfile)) redirec = paste(">", outfile)
# system(paste(tabixcmd, gzpath, chrom, recfilt, redirec), intern=TRUE)
#}
#setClass("metaVCF", contains="character")
#
#getMetaVCF = function (gzcon, maxnlines = 100, final = "^#CHROM") 
#{
#    if (!is(gzcon, "gzfile")) 
#        stop("gzcon must inherit from gzfile; see ?gzfile")
#    if (summary(gzcon)$opened != "opened") 
#        stop("please pass an opened connection")
#    meta = NULL
#    for (i in 1:maxnlines) {
#        meta = c(meta, tmp <- readLines(gzcon, n = 1))
#        if (length(grep(final, tmp)) > 0) 
#            break
#    }
#    close(gzcon)
#    new("metaVCF", meta)
#}
#
#setMethod("show", "metaVCF", function(object) cat(object, sep="\n"))
#
#
#sampleIDs = function(metavec, ndrop=9) 
#    strsplit( metavec[length(metavec)], "\t" )[[1]][-c(1:ndrop)]
#
##setMethod("sampleNames", c("metaVCF", "missing"), function(object, ...)
##  .sampleNames( object, ndrop=9 ) )
#
#setGeneric("vcf2sm", function(tbxfi, ..., gr, nmetacol) standardGeneric("vcf2sm"))
#setMethod("vcf2sm", c("TabixFile", "GRanges", "integer"), 
#   function( tbxfi, ..., gr, nmetacol=9) {
#     # get metadata, terminated by ^#CHROM
#     meta = NULL
#     while( length( tmp <- yieldTabix(tbxfi, yieldSize=1) ) > 0) {
#         meta = c(meta, tmp)
#         if (length(grep("^#CHROM", tmp))>0) break
#         }
#     sampids = sampleIDs(meta, ndrop=nmetacol)
#     Rsamtools:::close.TabixFile(tbxfi)
#     Rsamtools:::open.TabixFile(tbxfi)
#     chunk = scanTabix(tbxfi, param=gr)  # list of vectors of strings, one list elem per range in gr
#     out = list()
#     trk = 0
#     for (i in 1:length(chunk)) {
#       out[[i]] = lapply( chunk[[i]], function(x) { trk <<- trk+1; if (options()$verbose) if (trk %% 100 == 0) cat(trk);
#                    parseVCFrec( x, nmetacol=nmetacol, makelocpref="chr" ) })
#       }
#     out = unlist(out, recursive=FALSE)
#     rsid = sapply(out, "[[", "id")
#     nsnp = length(out)
#     mat = matrix(as.raw(0), nr=length(sampids), ncol=nsnp)
#     for (i in 1:nsnp) mat[,i] = out[[i]]$calls
#     rownames(mat) = sampids
#     colnames(mat) = rsid
#     Rsamtools:::close.TabixFile(tbxfi)
#     new("SnpMatrix", mat)
#})
#
#.vcf2sm = function(gzpath, chrom, tabixcmd = "tabix", nmetacol=9, verbose=FALSE,
#gran=10000, metamax=100, makelocpref="chr") {
# #require(snpMatrix)
# con = gzfile(gzpath, "r")
# mm = getMetaVCF( con, maxnlines=metamax )
# sampids = sampleIDs(mm, ndrop=nmetacol)
# on.exit({close(con); close(fpipe)})
# fpipe = filterVCF( gzpath, chrom, return.pipe = TRUE, tabixcmd = tabixcmd )
# out = list()
# i = 1
# while ( length(tmp <- readLines(fpipe, n=1))>0) {
#  out[[i]] = parseVCFrec( tmp, nmetacol=nmetacol, makelocpref=makelocpref )
#  if (verbose & (i%%gran)==0) cat(i)
#  i = i+1
#  }
# rsid = sapply(out, "[[", "id")
# nsnp = length(out)
# mat = matrix(as.raw(0), nr=length(sampids), ncol=nsnp)
# for (i in 1:nsnp) mat[,i] = out[[i]]$calls
# rownames(mat) = sampids
# colnames(mat) = rsid
# new("SnpMatrix", mat)
#}
#
#vcf2smTXT = function (txtpath, meta, nmetacol = 9, verbose = FALSE, gran=10000) 
#{
#    #require(snpMatrix)
#    mm = meta
#    sampids = sampleIDs(mm, ndrop = nmetacol)
#    on.exit(close(fpipe))
#    fpipe = pipe(paste("cat", txtpath), open="r")#  filterVCF(gzpath, chrom, return.pipe = TRUE, tabixcmd = tabixcmd)
#    out = list()
#    i = 1
#    while (length(tmp <- readLines(fpipe, n = 1)) > 0) {
#        out[[i]] = parseVCFrec(tmp, nmetacol = nmetacol)
#        if (verbose & (i%%gran) == 0)
#            cat(i)
#        i = i + 1
#    }
#    rsid = sapply(out, "[[", "id")
#    nsnp = length(out)
#    mat = matrix(as.raw(0), nr = length(sampids), ncol = nsnp)
#    for (i in 1:nsnp) mat[, i] = out[[i]]$calls
#    rownames(mat) = sampids
#    colnames(mat) = rsid
#    new("snp.matrix", mat)
#}
#
##vcf2metaloc = function(gzpath, chrom, tabixcmd = "tabix", nmetacol=9, verbose=FALSE,
##gran=10000) {
## #require(snpMatrix)
## mm = getMetaVCF( gzfile(gzpath, "r") )
## sampids = sampleIDs(mm, ndrop=nmetacol)
## on.exit(close(fpipe))
## fpipe = filterVCF( gzpath, chrom, return.pipe = TRUE, tabixcmd = tabixcmd )
## out = list()
## i = 1
## while ( length(tmp <- readLines(fpipe, n=1))>0) {
##  out[[i]] = parseVCFrec( tmp, nmetacol=nmetacol )
##  if (verbose & (i%%gran)==0) cat(i)
##  i = i+1
##  }
## rsid = sapply(out, "[[", "id")
## locs = as.numeric(sapply(out, "[[", "loc"))
## ref = sapply(out, "[[", "ref")
## alt = sapply(out, "[[", "alt")
## depth = sapply(out, "[[", "depth")
## depth = gsub(".*DP=", "", depth)
## depth = as.numeric(depth)
## require(GenomicRanges)
## rng = GRanges(seqnames=paste("chr", chrom, sep=""), IRanges(start=locs, width=1))
## names(rng) = rsid
## elementMetadata(rng) = DataFrame(ref=ref, alt=alt, depth=depth)
## rng
##}
#
parseVCFrec = function(rec, nmetacol=9, makelocpref="chr" ) {
 vec = strsplit(rec, "\t")[[1]]
 meta = vec[1:nmetacol]
 calls = vec[-c(1:nmetacol)]
 nalt = strsplit(calls, "")
 nums = lapply(nalt, "[", c(1,3))  # extract the call components
 hasmiss = which(sapply(nums, function(x) any(x == ".")))
 nalt = sapply(nums, function(x) 2-sum(x=="0"))  # this is correct only for diallelic locus; note in doc
 if (length(hasmiss)>0) nalt[hasmiss] = -1
 nalt = nalt+1
 chr = meta[1]
 id = meta[3]
 loc = meta[2]
 if (id == "." ) id = paste(makelocpref, chr, ":", loc, sep="")
 list(chr=chr, id=id, loc=loc, ref=meta[4], alt=meta[5], depth=meta[8],
   calls=as.raw(nalt))
}

setGeneric("vcf2sm", function(tbxfi, ..., gr, nmetacol) standardGeneric("vcf2sm"))
setMethod("vcf2sm", c("TabixFile", "GRanges", "integer"), 
   function( tbxfi, ..., gr, nmetacol=9) {
     # get metadata, terminated by ^#CHROM
     # get fresh connection
     if (Rsamtools:::isOpen(tbxfi)) Rsamtools:::close.TabixFile(tbxfi)
     Rsamtools:::open.TabixFile(tbxfi)
     head = scanVcfHeader(tbxfi)  # you are positioned at data
     sampids = head@samples # head[[1]][["Sample"]] # some reflectance
     chunk = scanTabix(tbxfi, param=gr)  # list of vectors of strings, one list elem per range in gr
     out = list()
     trk = 0
     for (i in 1:length(chunk)) {
       if (length(chunk[[i]]) == 0) next
       out[[i]] = lapply( chunk[[i]], function(x) { trk <<- trk+1; if (options()$verbose) if (trk %% 100 == 0) cat(trk);
                    parseVCFrec( x, nmetacol=nmetacol, makelocpref="chr" ) })
       }
     out = unlist(out, recursive=FALSE)
     if (length(out) == 0) return(NULL)
     rsid = sapply(out, "[[", "id")
     nsnp = length(out)
     mat = matrix(as.raw(0), nr=length(sampids), ncol=nsnp)
     for (i in 1:nsnp) mat[,i] = out[[i]]$calls
     rownames(mat) = sampids
     colnames(mat) = rsid
     Rsamtools:::close.TabixFile(tbxfi)
     new("SnpMatrix", mat)
})

