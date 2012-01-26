tname = function() gsub(".*/", "", tempfile())

cisTo = function (sr, gr, rad = 2e+06) 
{
    ac = function(x) as(x, "character")
    if (length(unique(ac(seqnames(sr)))) > 1) stop("need sr to represent only one chromosome")
    if (!all(unique(ac(seqnames(sr))) ==  unique(ac(seqnames(gr))))) {warning("sr and gr represent different chromosomes"); return(NULL)}
    rope = sr + rad
    fo = findOverlaps(rope, gr)
    ss = split(fo@matchMatrix[,2], fo@matchMatrix[,1])
    uu = unique(fo@matchMatrix[,1])
    nsr = names(sr)[uu]
    gn = names(gr)
    ans = lapply(ss, function(x) gn[x])
    names(ans) = nsr
    ans
}

getTransScoresOneSNP = function(snpname, mgr, cisToList=NULL, ffind=1, n=200, searchfac=4 ) {
#
# one SNP, same-chromosome application
# cisToList has one element per SNP, names of probes within a given radius (use "cisTo" to generate)
# ffind is the manager component harboring SNP, n is number of scores to retrieve, searchfac is
# factor by which n is multiplied to obtain scores prior to eliminating cis scores
#
 if (length(snpname)>1) stop(paste("supply only a single SNP name, got", length(snpname)))
 sn = rsid( snpname )
 if (!is.null(cisToList) && !(snpname %in% names(cisToList))) stop("snpname not in cisToList")
 t1 = topFeats(sn, mgr=mgr, useSym=FALSE, ffind=ffind, n=n*searchfac )
 drops = NULL
 if(!is.null(cisToList)) drops = intersect(cisToList[[sn]], names(t1))
 if (length(drops)==0) return(sort(t1, decreasing=TRUE)[1:n])
 t1 = t1[ -match(drops, names(t1)) ]
 ot1 = order(t1, decreasing=TRUE)[1:min(c(length(t1),n))]
 sco1 = t1[ot1]
 n1 = names(t1)[ot1]
 list(scores=sco1, genes=n1)
}

#getTransScores = function(mgr, cisToList=NULL, ffind=1, 

#if (!exists("t1")) load("t1.rda")
#gg = getTransScores( rsid("rs4814683"), mgr=t1, cisToList=fff )
