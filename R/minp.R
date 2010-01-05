minp = function(path, BUFSIZE=100000, breakat=Inf, freqdump=10000,
  filefun=gzfile) {
   owarn = options()$warn
   ans = rep(NA, BUFSIZE)
   minind = rep(NA, BUFSIZE)
   rsid = rep(NA, BUFSIZE)
   on.exit(close(fi))
   nl = 0; fi = filefun(path, "r");
   hline = NULL
   while( !inherits(try(tmp <- readLines(fi, n=1, ok=FALSE)), "try-error") ) {
     if (nl == 0 & is.null(hline)) {  # get header record and don't increment nl
        hline = strsplit(tmp, " ")[[1]]
        next
        }
     if (nl %% freqdump == 0 && options()$verbose) cat(nl) # monitor
     nl = nl+1  # bump line count
     if (nl %% BUFSIZE == 0) {       # add buffer if needed
            ans = c(ans, rep(NA, BUFSIZE))  
            minind = c(minind, rep(NA, BUFSIZE))  
            rsid = c(rsid, rep(NA, BUFSIZE)) 
            if (options()$verbose) cat("growing buffer...\n")
            }
     txt = strsplit(tmp, " ")[[1]]
     rsid[nl] = as.numeric(txt[2])
     options(warn=0)
     chisqs = as.numeric(strsplit(tmp, " ")[[1]][-c(1,2)])
     options(warn=owarn)
     if (all(is.na(chisqs))) {
         ans[nl] = NA
         next
         }
     if (nl >= breakat) break
     N = sum(!is.na(chisqs))
     rawp = pmin(1, 2*(1-pchisq(chisqs,1)))
     ans[nl] = min(rawp)
     minind[nl] = which.min(rawp)
# following chunk is intended to deal with multiplicity and correlation
# through Benjamini-Yekutieli algorithm -- but passing for now
   #  adjp.tmp = mt.rawp2adjp(rawp, "BY")
   #  adjp = adjp.tmp$adjp[,"BY"][order(adjp.tmp$index)]
   #  if (all(adjp >=1 )) {
   #    ans[nl] =NA
   #    next
   #    }
   #  fish = sum(-2*log(adjp),na.omit=TRUE) # Chisq w 2N DF
   #  print(fish)
   #  print(N)
   #  fishp = 1-pchisq(fish, 2*N)
   #  ans[nl] = abs(qnorm(fishp))
    }
    list(rsid=rsid[1:nl], ans=ans[1:nl], minind=minind[1:nl], genes=hline[-1])
}



locreport = function (winfo, chr) 
{
    if (!is.numeric(chr)) 
        stop("chr must be in 1:24")
    require(GGBase)
    locs = snpLocs.Hs(chrnum(chr), rsid(paste("rs", winfo[["rsid"]], 
        sep = "")))
    ldf = data.frame(reflocs = locs[2, ], rsid = paste("rs", 
        locs[1, ], sep = ""))
    gns = gsub("\"", "", winfo[["genes"]])
    picks = gns[winfo[["minind"]]]
    ac = as.character
    indf = data.frame(rsid = ac(paste("rs", winfo[["rsid"]], 
        sep = "")), mlogp = -log10(winfo[["ans"]]), mind = winfo[["minind"]], 
        bestgene = ac(picks), stringsAsFactors = FALSE)
    fulldf = merge(ldf, indf, by = "rsid", all.x = TRUE)
    fulldf$rsid = ac(fulldf$rsid)
    fulldf
}

wgtinfo2browser = function(winfo, chr, start, end, trname="newtrack",
  bsession=NULL, vlim=c(0,16),  ...) {
  require(rtracklayer, quietly=TRUE)
  lrep = locreport(winfo, chr)
  lrepc = na.omit(lrep)
  lrepco = lrepc[ order(lrepc$reflocs), ]
  space = paste("chr", chr, sep="")
  tsco = lrepco$mlogp
  tsco[!is.finite(tsco)] = max(tsco[is.finite(tsco)]) # winsorize at top
  rr = RangedData(IRanges(start=lrepco[,2], end=lrepco[,2]),
      score = tsco, space=space)
  if (is.null(bsession)) ee = browserSession("UCSC")
  else ee = bsession
  track(ee, trname, viewLimits=vlim) = rr
  browserView(ee, GenomicRanges(start=start, end=end, chrom=space),
    track=c("ruler", "cytoBand"), full=trname, ...)
}
 
top4 = function (x, sms) 
{
# assumes you have split a locreport output by gene
    par(mfrow = c(2, 2))
    for (i in 1:4) plot_EvG(probeId(names(x)), rsid(x[[1]][i, 
        "rsid"]), sms, main = paste("-log10p=", round(x[[1]][i, "mlogp"], 
        3)))
}

setClass("maxchisq", contains="list")
setMethod("show", "maxchisq", function(object) {
 cat("GGtools maxchisq structure.\n")
 cat("The call was:\n") 
 print(object$theCall)
 cat("The original call in multffManager was:\n")
 print(object$mgrcall)
 cat("Excerpt:\n")
 print(lapply(object[c("maxchisq", "bestgenes")], function(x) head(x[[1]])))
})

maxchisq = function(mgr, nchr=length(mgr$fflist), ncores=2) {
 theCall = match.call()
 mgrcall = mgr$call
 require(multicore)
 maxchisq = mclapply( mgr$fflist[1:nchr],  function(x) apply(x[], 1, max)/mgr$shortfac)
 bestgenes = mclapply( mgr$fflist[1:nchr],  function(x) colnames(x)[apply(x[], 1, which.max)]) 
 new("maxchisq", list(maxchisq=maxchisq, df=mgr$df, bestgenes=bestgenes, theCall=theCall, mgrcall=mgrcall))
}

setGeneric("min_p_vals", function(mcs, mtcorr, type) standardGeneric("min_p_vals"))
setMethod("min_p_vals", c("maxchisq", "character", "character"), function(mcs, mtcorr, type) {
 pv = lapply(mcs$maxchisq, function(x) pmin(1, 2*(1-pchisq(x, mcs$df))))
 mtcorrp = function(x, mtcorr) {
   tmp = mt.rawp2adjp(x, mtcorr)
   tmp$adjp[ order(tmp$index), mtcorr ]
 }
 if (mtcorr == "none") adjpv = pv
 else {
   require(multtest)
   if (type == "chr_specific")
     adjpv = lapply( pv, function(x) mtcorrp(x, mtcorr))
   else if (type=="global") {
     ulp = unlist(pv)
     uln = unlist(mcs$bestgenes)
     adjpv = mtcorrp(ulp, mtcorr)
     names(adjpv) = uln
#     return(adjpv) # returns one named vector
     anslist = list()
     for (i in 1:length(mcs$bestgenes)) {
       anslist[[i]] = adjpv[ mcs$bestgenes[[i]] ] # restore chromosomal list structure
     }
     names(anslist) = names(mcs$bestgenes)
     return(anslist)
   }  
 }
 for (i in 1:length(adjpv))
   names(adjpv[[i]]) = mcs$bestgenes[[i]]
 adjpv
})



