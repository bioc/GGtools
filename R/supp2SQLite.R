supp2SQLite = function (df, tabname, fn, ...) 
{
    require(RSQLite)
    require(Biostrings)
    require(org.Hs.eg.db)
    imap = names(IUPAC_CODE_MAP)
    names(imap) = as.character(IUPAC_CODE_MAP)
    aas = as.character(df$Assignment)
    aas = gsub("\\.|/", "", aas)
    aas[nchar(aas) != 2] = "ACGT"
    aasm = imap[aas]
    rsnn = as.numeric(gsub("rs", "", rownames(df)))
    cnum = gsub("chr", "", df$Chromosome)
    if (length(unique(cnum)) != 24) 
        stop("number of chromosome tokens must equal 24")
    cnum[cnum == "X"] = 23
    cnum[cnum == "Y"] = 24
    cnum = as.numeric(cnum)
    locs = df$Position
    slocs = split(locs, cnum)
#    mxs = sapply(slocs, max)
#    cmxs = cumsum(as.numeric(mxs))
    cmxs = cumsum(as.numeric(org.Hs.egCHRLENGTHS[1:24]))
    for (i in 2:24) slocs[[i]] = slocs[[i]] + cmxs[i - 1]
    cumloc = unlist(slocs)
    newdf = data.frame(rsid = as.integer(rsnn), alleles = as.character(aasm), 
        chrnum = as.integer(cnum), loc = cumloc)
    newdf$alleles = as.character(newdf$alleles)
    #df2nc(newdf, fn, ...)
    drv = dbDriver("SQLite")
    con = dbConnect(drv, fn)
    cat("writing table ...\n")
    dbWriteTable(con, tabname, newdf, row.names=FALSE)
    cat("done.")
    invisible(newdf)
}
