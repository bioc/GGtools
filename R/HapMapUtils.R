HMworkflow = function(gzfn,emat,pd,mi,anno) {
# assumes gzfn is pathname of a gzipped HapMap file
# emat is rownames/columnames matrix of exprs
# pd is a phenoData structure
 savo = options()
 on.exit(options(savo))
 options(verbose=TRUE)
 require(GGtools)
 rac = thinHM2rac(gzfn)
 racm = rac$rarecount[, colnames(emat)] 
 new("racExSet", exprs=emat, racs=racm, rarebase=rac$rarebase,
   SNPalleles=rac$alleles, phenoData=pd, experimentData=mi,
   annotation=anno)
}


HM2rac = function (fn, comment.char = "", kppref = "NA") 
{
# SLOW!  see thinHM2rac below
    df = read.table(fn, comment.char = comment.char, h = TRUE)
    allel = as.character(df[, "SNPalleles"])
    kpcol = grep(kppref, names(df))
    df = df[, c(1, kpcol)]
    names(df)[1] = "rsnum"
    BAG = lapply(df, as.character)
    nsnp = nrow(df)
    nsam = ncol(df) - 1
    cmat = matrix(" ", nrow = nsnp, nc = nsam)
    for (i in 1:nsam) cmat[, i] = BAG[[i + 1]]
    names(allel) = rownames(cmat) = BAG[[1]]
    colnames(cmat) = names(df)[-1]
    rr = apply(cmat, 1, countRare)
    rrr = apply(cmat, 1, getRare)
    vdf = data.frame("sampID"=colnames(cmat))
    vmdf = data.frame(labelDescription="sample ID in CEPH system")
    class(vmdf[,1]) = "character"
    vmdf[1,1] = "sample ID in CEPH system"
    rownames(vmdf) = "sampID"
    adf = new("AnnotatedDataFrame", data=vdf, varMetadata=vmdf)
    list(raremat = t(rr), alleles=allel, rareallele = rrr, anno=adf)
}


thinHM2rac = function(gzfn) {
#
# here we scan in the gzipped data line by line and
# build up the matrix
#
 stats = system(paste("gunzip -c", gzfn, "|wc"), intern=TRUE)
 stats = as.numeric(strsplit(stats, "\ +")[[1]])
 ntokpline = stats[3]/stats[2]
 nline = stats[2]-1 # header
 nind = ntokpline - 11
 rarecount = matrix(NA, nr=nline, nc=nind)
 alleles = rep(NA, nline)
 rarebase = rep(NA, nline)
 rsnum = rep(NA, nline)
 ff = gzfile(gzfn)
 open(ff, "r")
 hd = scan(ff, "", n=ntokpline, quiet=TRUE)
 snames = hd[-c(1:11)]
 cat(paste(nline, "lines to process\n"))
 for (i in 1:nline)
   {
   if (i %% 500 == 0) cat(i)
   tmp = scan(ff, "", n=ntokpline, quiet=TRUE)
   rarebase[i] = getRare(tmp[-c(1:11)])
   rarecount[i,] = countRare(tmp[-c(1:11)])
   rsnum[i] = tmp[1]
   alleles[i] = tmp[2]
   }
 rownames(rarecount) = rsnum
 colnames(rarecount) =  snames
 list(rarecount=rarecount, rarebase=rarebase, 
       alleles=alleles)
}
 
