
cleanup = function (x, tag = "snp.", convert=FALSE) 
{
    tmp = gsub("\\\"", "", gsub(tag, "", x))
    if (!convert) return(tmp)
    as.numeric(tmp)
}

sgf.unsafe = function (x) 
{
#
# output of tabixScan of cisRun gff3 is converted to data frame
#
    if (length(x) == 0) return(list(snp=NULL, gene=NULL, fdr=NULL))
    z = strsplit(x, "\t")
    chr = sapply(z, "[", 1)
    sco = sapply(z, "[", 6)
    z9 = sapply(z, "[", 9)
    zz = strsplit(z9, "; ") # always trailing space?  or do in cleanup
    snp = sapply(zz, "[", 1)
    loc = sapply(zz, "[", 2)
    gene = sapply(zz, "[", 6)
    fdr = sapply(zz, "[", 5)
    data.frame(snp=cleanup(snp), gene=cleanup(gene, "probeid."), 
       fdr=cleanup(fdr, "fdr.", convert=TRUE), score=as.numeric(sco),
       chr=chr, loc=cleanup(loc, "snplocs.", convert=TRUE), 
       stringsAsFactors=FALSE, row.names=NULL)
}



tscan2df = function( tscan ) {
 require(foreach)
 #library(doParallel)
 #registerDoParallel(cores=ncores)
 lk = foreach(i=1:length(tscan)) %dopar% sgf(tscan[[i]])
 names(lk) = NULL
 do.call(rbind, lk)
}



ceuCon19 = function() {
 fi = system.file("tabix/ceu1kgv_npc19.gff3.gz", package="piFDR2")
 TabixFile(fi)
}

unquote = function(tag) gsub("\\\"", "", tag)

sgf.slow = function (x) 
{
#
# output of tabixScan of cisRun gff3 is converted to data frame
#
    if (length(x) == 0) return(list(snp=NULL, gene=NULL, fdr=NULL))
    z = strsplit(x, "\t")
    chr = sapply(z, "[", 1)
    sco = sapply(z, "[", 6)
    z9 = sapply(z, "[", 9)
#
# code below is not safe if NAs are not propagated to GFF3
#
#    zz = strsplit(z9, "; ") # always trailing space?  or do in cleanup
#    snp = sapply(zz, "[", 1)
#    loc = sapply(zz, "[", 2)
#    gene = sapply(zz, "[", 6)
#    fdr = sapply(zz, "[", 5)
    atts2get = rep(" ", 4)
    names(atts2get) = c("snp", "loc", "gene", "fdr")
  zz = strsplit(z9, "; ") # always trailing space?  or do in cleanup

    aroundeq = lapply(zz, strsplit, "=")
    stuff = lapply(aroundeq, sapply, "[", 2)
    nstuff = lapply(aroundeq, sapply, "[", 1)
    stuff = lapply(1:length(stuff), 
         function(x) {
           names(stuff[[x]])=nstuff[[x]]
           stuff[[x]]
         })
    stuff = lapply(stuff, unquote)
    

    data.frame(snp=sapply(stuff, "[", "snp"), gene=sapply(stuff, "[", "probeid"),
       fdr = as.numeric(sapply(stuff, "[", "fdr")),
       score=as.numeric(sco),
       chr=chr, 
       loc = as.numeric(sapply(stuff, "[", "snplocs")),
       stringsAsFactors=FALSE, row.names=NULL)
}

getsfp = function(x) x[grep("snp=|snplocs|fdr|probeid", x)]  # guaranteed to be present unless NA

sgf = function (x)
{
    if (length(x) == 0)
        return(list(snp = NULL, gene = NULL, fdr = NULL))
    z = strsplit(x, "\t")
    chr = sapply(z, "[", 1)
    sco = sapply(z, "[", 6)
    z9 = sapply(z, "[", 9)
    zz = strsplit(z9, "; ")
    zzpull = sapply(zz, getsfp)
    stopifnot(is.matrix(zzpull))  # means some of snp=, fdr, probeid not present in some record
    snp = unquote(sapply(strsplit(zzpull[1,], "="), "[", 2))
    loc = as.numeric(sapply(strsplit(zzpull[2,], "="), "[", 2))
    fdr = as.numeric(sapply(strsplit(zzpull[3,], "="), "[", 2))
    fdrtag = sapply(strsplit(zzpull[3,], "="), "[", 1)
    stopifnot(all(fdrtag=="fdr"))  # sanity check
    gene = unquote(sapply(strsplit(zzpull[4,], "="), "[", 2))
    data.frame(snp = snp, gene = gene,
         fdr = fdr,
         score = as.numeric(sco), chr = chr, loc = loc,
         stringsAsFactors = FALSE, row.names = NULL)
}

#[[2]]
#[1] "chr17"                                                                                                                                                                                                
#[2] "rtracklayer"                                                                                                                                                                                          
#[3] "sequence_feature"                                                                                                                                                                                     
#[4] "828"                                                                                                                                                                                                  
#[5] "828"                                                                                                                                                                                                  
#[6] "0.04"                                                                                                                                                                                                 
#[7] "-"                                                                                                                                                                                                    
#[8] "."                                                                                                                                                                                                    
#[9] "snp=\"rs62053745\"; snplocs=828; ests=0; se=0.01; fdr=0.98530721059859; probeid=\"EnnUlQgJcHk1dN19nk\"; MAF=0.151898734177215; dist.mid=17887.5; permScore_1=1.97; permScore_2=1.42; permScore_3=2.63"

getall = function(x) x[grep("snp=|snplocs|fdr|probeid|MAF|permScore_1|permScore_2|permScore_3", x)]  # guaranteed to be present unless NA
# ests and se not guaranteed

fullparse = function(x) {
    if (length(x) == 0)
        return(GRanges())
    z = strsplit(x, "\t")
    chr = sapply(z, "[", 1)
    sco = sapply(z, "[", 6)
    snpstart = sapply(z, "[", 4)
    z9 = sapply(z, "[", 9)
    zz = strsplit(z9, "; ")
    zzpull = sapply(zz, getall)
    stopifnot(is.matrix(zzpull))  # means some of snp=, fdr, probeid not present in some record
    snp = unquote(sapply(strsplit(zzpull[1,], "="), "[", 2))
    fdr = as.numeric(sapply(strsplit(zzpull[3,], "="), "[", 2))
    fdrtag = sapply(strsplit(zzpull[3,], "="), "[", 1)
    stopifnot(all(fdrtag=="fdr"))  # sanity check
    gene = unquote(sapply(strsplit(zzpull[4,], "="), "[", 2))
    maf = as.numeric(sapply(strsplit(zzpull[5,], "="), "[", 2))
    maftag = sapply(strsplit(zzpull[5,], "="), "[", 1)
    perm1 = as.numeric(sapply(strsplit(zzpull[6,], "="), "[", 2))
    perm2 = as.numeric(sapply(strsplit(zzpull[7,], "="), "[", 2))
    perm3 = as.numeric(sapply(strsplit(zzpull[8,], "="), "[", 2))
    stopifnot(all(maftag=="MAF"))  # sanity check
    GRanges(seqnames=chr, IRanges(as.numeric(snpstart),width=1), snp = snp, gene = gene,
         fdr = fdr,
         score = as.numeric(sco), 
         maf = maf, perm1=perm1, perm2=perm2, perm3=perm3)
}

tscan2gr = function( tscan ) {
 require(foreach)
 #library(doParallel)
 #registerDoParallel(cores=ncores)
 lk = foreach(i=1:length(tscan), .combine=c) %dopar% fullparse(tscan[[i]])
 names(lk) = NULL
 lk
}

genemodel = function(sym) {
 require(Homo.sapiens)
 egid = select(Homo.sapiens, 
     keytype="SYMBOL",
     keys=sym, 
     columns=c("GENEID", "CHR", "TXID", "TXSTART", "TXEND", "TXSTRAND",
       "EXONID", "EXONNAME", "EXONSTART", "EXONEND"))
 GRanges(seqnames=paste0("chr", egid$CHR[1]), strand=egid$TXSTRAND,
   IRanges(egid$EXONSTART, egid$EXONEND), symbol=egid$SYMBOL,
      entrezid=egid$GENEID, transcript=egid$TXID, exon=egid$EXONID, chromosome=paste0("chr", egid$CHR[1]))
}

#
#
# orm = genemodel("ORMDL3")
# t1 autoplot(orm, layout="linear")

queryCrg = function (attr, gr, con )
{
#
# obtain attributes named in attr after scanTabix(con, param=gr)
# scan may be conducted in parallel using foreach/doParallel
#
    dat = scanTabix(con, param = gr)
    df = tscan2df(dat)
    df[, attr]
}

crg2gr = function (gr, con )
{
#
# obtain attributes named in attr after scanTabix(con, param=gr)
# scan may be conducted in parallel using foreach/doParallel
#
    dat = scanTabix(con, param = gr)
    tscan2gr(dat)
}
scoresCis = function( sym="ORMDL3", cisRun, cisannopk="lumiHumanAll.db", 
    radius=100000,
    pad=1000, 
    txScore = function(x) -log10(x+(1e-5)), ylim=c(0,4), genometag="hg19", plot.it=TRUE,
    laxistag = "-log10 FDR: ", ...) {
#
# on the basis of a data.table output of ciseqByCluster, display transformed scores
#
 stopifnot(is(cisRun, "data.table"))
 m = genemodel(sym)
 require(cisannopk, character.only=TRUE)
 pid = get(sym, revmap(get(paste0(gsub(".db", "", cisannopk), "SYMBOL"))))
 matcher = cisRun[ which(cisRun$probeid %in% pid), list(seqnames,fdr,start) ] 
 runnear = GRanges( matcher$seqnames[1], IRanges(matcher$start, width=1) )
 sc = txScore( as.numeric(matcher$fdr) )
 stopifnot(length(sc)==length(runnear))
 G = DataTrack(runnear, ylim=ylim, name=paste0(laxistag, sym), genome=genometag)
 values(G) = matrix(sc, nrow=1)
 ans = list(dtrack=G,
     grt=GeneRegionTrack(m, showId=TRUE), 
     gat=GenomeAxisTrack(name=seqlevels(m), showTitle=TRUE, col.title="gray"))
 if (plot.it) plotTracks(ans, ...)
 invisible(ans)
}
