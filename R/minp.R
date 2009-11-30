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

