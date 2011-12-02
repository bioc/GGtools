eqtlEstimates = function (smlSet, rhs = ~1 - 1, runname = "fooe", targdir = "fooe", 
    geneApply = lapply, chromApply = lapply, shortfac = 100, 
    checkValid = TRUE, saveSummaries = TRUE, uncert = TRUE, family, 
    genegran = 50, prefilter = dropMonomorphies, ...) 
{
    theCall = match.call()
    if (is.function(prefilter)) 
        smlSet = prefilter(smlSet)
    if (checkValid) {
        tmp = validObject(smlSet)
    }
    if (missing(family)) 
        family = "gaussian"
    geneindex <- 1
    sess = sessionInfo()
    fnhead = paste(targdir, "/", runname, "_", sep = "")
    geneNames = featureNames(smlSet)
    chrNames = names(smList(smlSet))
    ngenes = length(geneNames)
    nchr = length(chrNames)
    if (!file.exists(targdir)) 
        dir.create(targdir)
    summfflist = list()
    if (saveSummaries) {
        sumfn = paste(fnhead, chrNames, "_summ.ff", sep = "")
        if ("multicore" %in% search()) {
            summfflist = multicore::mclapply(1:length(chrNames), 
                function(i) ffSnpSummary(smList(smlSet)[[i]], 
                  sumfn[i], fac = shortfac))
        }
        else {
            for (i in 1:length(sumfn)) summfflist[[chrNames[i]]] = ffSnpSummary(smList(smlSet)[[i]], 
                sumfn[i])
        }
    }
    cres = chromApply(chrNames, function(chr) {
        snpdata = smList(smlSet)[[chr]]
        targff = paste(fnhead, "chr", chr, ".ff", sep = "")
        snpnames = colnames(snpdata)
        nsnps = ncol(snpdata)
        if (file.exists(targff)) 
            cat("attempting to overwrite ", targff, "...")
        store = ff(initdata = 0, dim = c(nsnps, ngenes, 2), dimnames = list(snpnames, 
            geneNames, c("est.", "s.e.")), vmode = "short", filename = targff)
        geneApply(geneNames, function(gene) {
            if (options()$verbose & geneindex%%genegran == 0) 
                cat(gene, "..")
            geneindex <- geneindex + 1
            if (options()$verbose & geneindex%%8 * genegran == 
                0) 
                cat("\n")
            ex = exprs(smlSet)[gene, ]
            fmla = formula(paste("ex", paste(as.character(rhs), 
                collapse = ""), collapse = " "))
            numans.full = snp.rhs.estimates(fmla, snp.data = snpdata, 
                data = pData(smlSet), family = family, uncertain = uncert, 
                ...)
            numans = sapply(numans.full, "[[", "beta")
            if (any(badests <- sapply(numans, is.null))) {
                numans[badests] = NA
                numans = unlist(numans)
		}
            numans.se = sapply(numans.full, "[[", "Var.beta")
            if (any(badses <- sapply(numans.se, is.null))) {
                numans.se[badses] = NA
                numans.se = sqrt(unlist(numans.se))
            }
            store[, gene, 1, add = TRUE] = shortfac * numans
            store[, gene, 2, add = TRUE] = shortfac * numans.se
            NULL
        })
        close(store)
        store
    })
    names(cres) = chrNames
    exdate = date()
    new("eqtlEstimatesManager", fflist = cres, call = theCall, 
        sess = sess, exdate = exdate, shortfac = shortfac, geneanno = annotation(smlSet), 
        df = 1, summaryList = summfflist)
}


setMethod("[", "eqtlEstimatesManager", function(x, i, j, k, ..., drop=FALSE) {
#
# ultimately this may not be exposed, serving only for deep
# testing, because a director database may be required for every
# manager
#

 if (!(k %in% c(1L,2L))) stop("3rd index must be 1L (for estimates) or 2L (for s.e.s)")
 m1 = snpIdMap( as(i, "character"), x )
#
# do not rebind i here
#
 ans = lapply(1:length(m1), function(mapi) fflist(x)[[names(m1)[mapi]]][ m1[[mapi]], 
    as(j, "character"), k, drop=FALSE]/shortfac(x))
 names(ans) = names(m1)
 ans
})

setMethod("[", c("eqtlEstimatesManager"), 
 function(x, i, j, k,..., drop=FALSE) {
 if (!missing(i) & !missing(j) & !missing(k)) {
   if (!is(i, "rsid")) stop("index i must be rsid instance")
   if (!is(j, "probeId")) stop("index j must be probeId instance")
   if (!is(k, "integer")) stop("index k must be integer")
     m1 = snpIdMap( as(i, "character"), x )
     ans = lapply(1:length(m1), function(mind) fflist(x)[[names(m1)[mind]]][ m1[[mind]],
        as(j, "character"), k, drop=FALSE]/shortfac(x))
   }
 else if (missing(i) & !missing(j) & !missing(k)) {
   if (!is(j, "probeId")) stop("index j must be probeId instance")
   if (!is(k, "integer")) stop("index k must be integer")
     ans = lapply(1:length(fflist(x)), function(mind) fflist(x)[[mind]][ ,
        as(j, "character"), k, drop=FALSE]/shortfac(x))
   }
 else stop("one of i (rsid instance), j (probeId instance) must be present along with k in (1L, 2L)")
 names(ans) = names(fflist(x))
 ans
 })

