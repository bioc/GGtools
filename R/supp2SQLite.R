supp2SQLite = function (df, tabname, fn, ...) 
{
    require(RSQLite)
    require(Biostrings)
    require(org.Hs.eg.db)
    indf = df
# first, order the rows of df ensuring that chromosomes are in order
#  a) get chrom tokens into numeric share
    cnum = gsub("chr", "", df$Chromosome)
    if (length(unique(cnum)) != 24) 
        stop("number of chromosome tokens must equal 24")
    cnum[cnum == "X"] = 23
    cnum[cnum == "Y"] = 24
# b) obtain order
    cnum = as.numeric(cnum)
    df  = df[order(cnum),]
    cnum = cnum[order(cnum)]
# map ambiguities
    imap = names(IUPAC_CODE_MAP)
    names(imap) = as.character(IUPAC_CODE_MAP)
    aas = as.character(df$Assignment)
    aas = gsub("\\.|/", "", aas)
    aas[nchar(aas) != 2] = "ACGT"
    aasm = imap[aas]
# get numeric representation of snpnum
    rsnn = as.numeric(gsub("rs|ss", "", rownames(df)))
    locs = df$Position
    slocs = split(locs, cnum)
# save offsets per chromosome
    offs = as.integer(sapply(slocs, min))
    cmxs = cumsum(as.numeric(org.Hs.egCHRLENGTHS[1:24]))
    for (i in 2:24) slocs[[i]] = slocs[[i]] + cmxs[i - 1]
    cumloc = unlist(slocs)
    noffs = rep(NA, length(cumloc))
    noffs[1:24] = offs
#    ord = order(as.integer(cnum), locs)
    newdf = data.frame(rsid = as.integer(rsnn), alleles = as.character(aasm), 
        chrnum = as.integer(cnum), loc = cumloc, offsets=noffs)
    newdf$alleles = as.character(newdf$alleles)
#    newdf = newdf[ord,]
#    newdf = cbind(newdf, offsets=as.integer(noffs))
    #df2nc(newdf, fn, ...)
    drv = dbDriver("SQLite")
    con = dbConnect(drv, fn)
    cat("writing table ...\n")
    dbWriteTable(con, tabname, newdf, row.names=FALSE)
    cat("done.")
    invisible(list(origdf=indf, sorteddf=df, newdf=newdf))
}
