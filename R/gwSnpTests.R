

setMethod("gwSnpTests", c("formula", "smlSet"),
  function( sym, sms, ...) {
    theCall = match.call()
    infmla = sym
#
# decode formula -- note that for a gene set we are essentially recursing
# so it is hard to factor this part out
#
    respObj = eval(sym[[2]]) # we know sym is a formula, sym[[2]] is dep var
    if (is(respObj, "phenoVar")) { pid = as.character(respObj@.Data) }
    else if (is(respObj, "genesym")) {
      annpack = sms@annotation
      require(annpack, character.only=TRUE)
      rmap = revmap( get(paste(gsub(".db", "", annpack), "SYMBOL", sep="")) )
      pid = AnnotationDbi::get( as(respObj, "character"), rmap )
      if (length(pid) == 0) stop(paste("cannot map", respObj, "in", annpack, sep=""))
      pid = intersect(pid, featureNames(sms))
      if (length(pid) > 1) {
        warning(paste("several probes/sets map to", respObj, "; using", pid[1], sep=""))
        print(pid) 
        pid = pid[1]
        }
      }
    else if (is(respObj, "probeId")) pid = respObj
    else stop("response in formula must be of class phenoVar, genesym, probeId, or GeneSet")
#
# at this point we have the featureName that we need
#
    pname = as.character(respObj)
    sym[[2]] = as.name(pname)  # replace the dependent variable spec in fmla
    if (!is(respObj, "phenoVar")) {
      assign(pname, exprs(sms)[pid,]) # expression phenotype genename
      alld = data.frame(get(pname), pData(sms))
      names(alld)[1] = pname
      }
    else alld = pData(sms)
    allsst = lapply( smList(sms), function(x) snp.rhs.tests(sym, family="gaussian",
        snp.data=x, data=alld, uncertain=TRUE))
    testType = "Gaussian"
## as of 8 july, we have data frames instead of snp.tests.single objects
## need to coerce
#    mksts = function(x) { 
#      new("snp.tests.single", chisq=cbind(`1 df`=x$Chi.squared, `2 df`=NA),
#          snp.names=rownames(x), N=x$Df.residual+x$Df, N.r2=numeric(0))
#    }
#    allsst = lapply(allsst, mksts)
#
#
    SI = new("SessionInfo", sessionInfo())
    new("gwSnpScreenResult", gene=respObj, psid=pid,
         annotation=sms@annotation, sessionInfo=SI, chrnum=names(smList(sms))[1],
         call=theCall, testType= testType, allsst)
    })

setMethod("show", "gwSnpScreenResult", function(object) {
 cat("gwSnpScreenResult for gene ", object@gene, " [probe ",
     object@psid, "]\n")
})

setMethod("topSnps", "gwSnpScreenResult", function(x, n=10) {
  ts.df = function (w, n = 10) {
   pp = p.value(w)
   sn = w@snp.names  # no accessor...  # don't need list access here
   opp = order(pp, decreasing=FALSE)
   spp = pp[ opp ]
   df = data.frame(p.val=spp)
   rownames(df) = sn[ opp ]
   df[1:n,,drop=FALSE]
   }
  lapply(x, ts.df, n=n)
})


genePosition = function (tok, genomeWide=FALSE, eglib="org.Hs.eg.db", annlib=NULL)
{
#
# internal only -- for deprecation -- i will not be jumping through these hoops
#
    tst = try(library(eglib, character.only=TRUE))
    egpref = gsub(".db", "", eglib)
    egSYM = get(paste(egpref, "SYMBOL", sep=""))
    egCHR = get(paste(egpref, "CHR", sep=""))
    egCHRLOC = get(paste(egpref, "CHRLOC", sep=""))
    egCHRLENGTHS = get(paste(egpref, "CHRLENGTHS", sep=""))
    if (!inherits(tst, "try-error")) {
        if (is(tok, "genesym")) 
            rmap = revmap(egSYM)
        else if (is(tok, "probeId")) {
            if (is.null(annlib)) stop("you must supply an annotation .db library name in annlib for a probeId token")
            require(annlib, character.only = TRUE, 
                quietly = TRUE)
            rmap = get(paste(gsub(".db", "", annlib), 
                "ENTREZID", sep = ""))
        }
        else {
            warning("tok is neither symbol nor probeID, we do not plot the location.")
            return(invisible(NULL))
        }
        egid = get(tok, rmap)
        ch = get(egid, egCHR)
        loc = get(egid, egCHRLOC)
        if (length(loc) > 1) {
            loc = loc[as.character(ch)][1]
        }
        if (length(loc) != 1) {
            warning("org.*.egCHRLOC has uninterpretable information for this gene; no tick at top of plot attempted.")
            return(invisible(NULL))
        }
        if (!genomeWide) return(abs(loc))
        chrl = egCHRLENGTHS
        chrbnd = cumsum(c(0,as.double(chrl[-length(chrl)])))
        if (ch == "X") 
            ch = 23
        else if (ch == "Y") 
            ch = 24
        else if (ch %in% c("Un", "MT")) 
            ch = 25
        else ch = as.numeric(ch)
        gpos = chrbnd[ch] + abs(loc)
        return(gpos)
    }
    warning("need org.*.eg.db for gene position.  no tick at top of plot attempted")
    return(invisible(NULL))
}

setMethod("plot", c("gwSnpScreenResult", "character"),  # y bound to location package
#
# for deprecation in 2012 -- use tracks!
#
  function(x, y=snplocsDefault(), noSmooth=FALSE, npts=500, ...) {
   allpv = snpStats:::p.value(x@.Data[[1]])
   rsn = x@.Data[[1]]@snp.names
   names(allpv) = rsn
   #if (is(x@.Data[[1]], "snp.tests.glm"))  # for new approach, snpMatrix2 > 1.1
   #      allpv = p.value(x@.Data[[1]])
   #else allpv = p.value(x@.Data[[1]]) #, y)
   kill = which(is.na(allpv))
   if (length(kill)>0) allpv = allpv[ -kill ]
   rsn = names(allpv)
   require(y, character.only=TRUE, quietly=TRUE)
   SNY = names(getSNPcount()) # new API -- shld suffice do.call(":::", list(y, "SEQNAMES"))
   if (!(x@chrnum %in% SNY))  x@chrnum = as.character(paste("chr", x@chrnum, sep=""))
   if (!(x@chrnum %in% SNY))  x@chrnum = as.character(gsub("chr", "ch", x@chrnum))
   if (!(x@chrnum %in% SNY))  stop("attempts to harmonize @chrnum of cwSnpScreenResult object with names(getSNPcount()) of snplocsDefault() failed")
#   loc = snpLocs.Hsapiens(rsn, x@chrnum, y) # may not match all

 require(y, character.only=TRUE)
 ldf = getSNPlocs(x@chrnum)
 if (nrow(ldf) == 0) stop("chrtok must be wrong")
 locs = ldf$loc
 names(locs) = paste("rs", as.character(ldf$RefSNP_id), sep="")
 chk = intersect(rsn, names(locs))
 if (length(chk) != length(rsid)) message(paste("NOTE: some SNP in rsid were not found in location db", y))
 if (length(chk) < 1) stop("no locations found for tested SNP")
 loc = locs[ chk ]

#   availRS = paste("rs", locstr["rsid",], sep="")
longnsubset = function (x, y)
{
    mm = match(y, names(x))
    x[mm]
}

   allpv = allpv[intersect(names(loc), names(allpv))]
   loc = loc[intersect(names(loc), names(allpv))]
#   allpv = longnsubset(allpv, availRS)
#   loc = locstr["loc",]
   if (noSmooth) plotf=plot
     else plotf=smoothScatter
   if (length(grep("resid", x@testType))>0) main = paste("resid", x@gene)
   else main=x@gene
   plotf(loc, -log10(allpv), main=main,
     xlab=paste("position on chr", x@chrnum),
     ylab=paste("-log10 p Gaussian LM [", 1, "df]", sep=""), pch=19, cex=.8, ...)
   #if (isCis(x))
        try(axis(3, at=genePosition(x@gene, annlib=x@annotation), col="red", lwd=2, label=" "))
})

