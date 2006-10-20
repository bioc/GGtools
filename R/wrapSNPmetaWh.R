wrapSNPmetaWh = function (rsn, chr, pos) 
{
    names(pos) = rsn
    spos = split(pos, chr)
    nch = names(spos)
    ok = which(!is.na(as.numeric(nch)))
    nok = max(as.numeric(nch[ok]))
    ospos = list()
    for (i in 1:nok) ospos[[i]] = spos[[as.character(i)]]
    names(ospos) = as.character(1:nok)
    if (length(spos) > nok) {
        n_nonnum = length(spos) - nok
        for (nn in nch[-ok]) ospos[[nn]] = spos[[nn]]
    }
    ends = cumsum(sapply(ospos, max))
    if (length(ospos) > 1) {
        for (i in 2:length(ospos)) ospos[[i]] = ospos[[i]] + 
            ends[i - 1]
    }
    whpos = unlist(ospos)
    mdf = data.frame(pos = whpos, strand = NA)
    rownames(mdf) = rsn
    mm = new.env()
    assign("meta", mdf, mm)
    new("snpMetaWhole", meta = mm, chrbounds = ends, chrlabs = names(ospos))
}
