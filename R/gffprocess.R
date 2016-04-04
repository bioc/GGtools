gffprocess = function(basename="fullyri100k", n_in=44, headpatt="_1A", tmpForSort="/freshdata/tmp") {
#
# given a collection of gff3 files representing chunked cis-eQTL search results, order them
# for concatenation and tabix indexing -- requires external classic unix grep and sort utilities
# along with tabix and bgzip
#
# we assume that -T argument to sort is followed by a folder path to be used for temporary files
# more modern sort utilities allow --temporary-directory= 
#
 allgff = dir(pattern="gff3$")
 stopifnot(length(allgff)==n_in)
 topind = grep(headpatt, allgff)
 stopifnot(length(topind)==1)
 tocat1 = allgff[-topind]
 tf1 = tempfile()
 tf2 = tempfile()
 run1 = system(paste0("cat ", paste(tocat1, collapse=" "), " | grep -v '^#' > ", tf1))
 stopifnot(run1 == 0)
 run2 = system(paste0("cat ", allgff[topind], " ", tf1, " > ", tf2 ))
 stopifnot(run2 == 0)
#!/bin/bash
 buildgrep = function(infile, outfile) {
#
# sorts and zips also
#
   templ = "(grep ^'#' %INF%; grep -v ^'#' %INF% | sort -T %TDIR% -k1,1 -k4,4n) | bgzip > %OUTF%"
   templ = gsub("%INF%", infile, templ)
   templ = gsub("%OUTF%", outfile, templ)
   templ = gsub("%TDIR%", tmpForSort, templ)
   templ
}
 outname = paste0(basename, ".gff3.gz")
 gcomm = buildgrep(tf2, outname) 
 run3 = system(gcomm) # now outname is bgzipped and sorted
 stopifnot(run3 == 0)
 run4 = system(paste0("tabix -p gff ", outname))
 stopifnot(run4 == 0)
 NULL
}
 
 
simpleTiling = function(ntile) {
 require(Homo.sapiens)
 hsi = seqinfo(Homo.sapiens)[paste0("chr", 1:22),]
 unlist(tileGenome(hsi, ntile=ntile))
 }

cgff2dt = function(gff3, tiling, addHitTest=TRUE, addcc878=TRUE ) {
#
# this program takes a genome-wide gff3 assembled by gffprocess after ciseqByCluster
# and a) assembles a data.table instance with the same content, b) computes the genome-wide FDR,
# c) builds a revised gff3 with gw
#
  requireNamespace("foreach")
  stopifnot(is(tiling, "GRanges"))
  basen = gsub(".gff3.gz", "", gff3)
  th = headerTabix(gff3)
  orderedChr = th$seqnames
  lgr = foreach(i=1:length(tiling)) %dopar% {
    gc()
    if (isTRUE(options()$gg.verbose)) cat(i)
    lk = try(import.gff3( gff3, which=tiling[i] ))
    if (inherits(lk, "try-error")) lk = NULL
    if (!is.null(lk)) lk = as.data.table(as.data.frame(lk))
    lk
  }
  bad = sapply(lgr, is.null)
  if (any(bad)) lgr = lgr[-which(bad)]
  ans = do.call(rbind, lgr)
  ans$snplocs = as.numeric(ans$snplocs)
  ans$ests = as.numeric(ans$ests)
  ans$se = as.numeric(ans$se)
  ans$oldfdr = as.numeric(ans$fdr)
  ans$MAF = as.numeric(ans$MAF)
  ans$dist.mid = as.numeric(ans$dist.mid)
  nperm = length(grep("permS", names(ans)))
  pnames = paste("permScore_", 1:nperm, sep="")
  for (i in 1:nperm)
    ans[[ pnames[i] ]] = as.numeric( ans[[ pnames[i] ]] )
  ans$mindist = as.numeric(ans$mindist)  # 
if (isTRUE(options()$gg.verbose)) print(date())

if (addcc878) {
 data(hmm878)
 eqr = GRanges(ans$seqnames, IRanges(ans$snploc, width=1))
 fo = findOverlaps(eqr, hmm878)
 nlev = length(unique(hmm878$name))
 chromcat878 = factor(rep("none", nrow(ans)), levels=c(unique(hmm878$name), "none"))
 chromcat878[ queryHits(fo) ] = factor(hmm878$name[subjectHits(fo)])
 ans$chromcat878 = chromcat878
}

if (addHitTest) {
 if (requireNamespace("gwascat")) {
   data(gwastagger)
   meqgr = GRanges(ans$seqnames, IRanges(ans$snploc, width=1))
   isgwashit = 1*(overlapsAny(meqgr, gwastagger) | ans$snp %in% gwastagger$tagid) # allow match by loc or name
   ans$isgwashit = isgwashit
   }
}

  ans$fdr = pifdr(ans$score, c(ans$permScore_1, ans$permScore_2, ans$permScore_3))
if (isTRUE(options()$gg.verbose)) print(date())
  obn = paste0(basen, "_dt")
  assign(obn, ans)
  save(list=obn, file=paste0(obn, ".rda"))
# prepare a file with the properly ordered genome-wide FDR
 gwfdr = ans$fdr
  gwch = ans$seqnames
  gwsnp = ans$snp
  gwprobeid = ans$probeid
#
# following code was used to verify synchronization to gff upon selecting by orderedChr
#
#  newdt = data.frame(gwch, gwsnp, gwprobeid, gwfdr, stringsAsFactors=FALSE)
#  newdt = split(newdt, gwch)
#  reviseFDR = file("newfdr.csv", open="a")
#  for (i in 1:length(orderedChr)) {
#      write.table(newdt[newdt$gwch == orderedChr[i],], file=reviseFDR,col.names=FALSE,row.names=FALSE,sep=",")
#  }
#  close(reviseFDR)
#
# use unix utilities to add the genome-wide FDR to the gff3
  sgwfdr = split(gwfdr, gwch)
  sgwfdr = unlist(sgwfdr[orderedChr])

    nn = tempfile()
    fdrt = file(nn, "w")

  writeLines(c(c(" ", " ", " "), paste("; gwfdr=", round(sgwfdr,6), sep="")), con=fdrt)
  close(fdrt)

  chkb = system(paste0("zcat ", gff3, " | paste -d '' - ", nn, " | bgzip > ", gsub(".gff3", "wmlf.gff3", gff3)))
  if(!(chkb == 0)) warning(paste("final system zcat-paste-bgzip returned non-zero: chkb=", chkb ))
  invisible(chkb)

}


