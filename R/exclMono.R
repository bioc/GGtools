
exclMono = function (res) 
{
# exclude monomorphic loci
    lu = apply(snps(res), 1, function(x) length(unique(x[!is.na(x)])))
    sn = snps(res)[lu > 1, ]
    make_racExSet(exprs(res), sn, rarebase(res), SNPalleles(res), 
        phenoData(res), experimentData(res), annotation(res))
}

