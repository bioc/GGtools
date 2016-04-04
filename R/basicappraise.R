
.discfmlas.demo = list(
  base=isgwashit~distcats+mafcats,
  basehmm=isgwashit~distcats+mafcats+chromcat878,
  allfac=isgwashit~distcats+mafcats+chromcat878+caddcats+fdrcats
)

.standardNames = c("seqnames", "start", "end", "width", "strand", "snp", 
 "snplocs", 
"score", "ests", "se", "fdr", "probeid", "MAF", "dist.mid", "mindist", 
"genestart", "geneend", "isgwashit", "chromcat878", "permScore_1", 
"permScore_2", "permScore_3", "PHRED")

appraise = function(dtab, discretize=TRUE, 
   reduceToSNP=TRUE, prefix, folder=paste0(prefix, "_APPROUT"), 
   discfmlas_in=GGtools:::.discfmlas.demo, txlist = list(
    distcats = function(x) {
     cut(x$mindist, c(-1, seq(0, 200001, 50000)))
     },
    fdrcats = function(x) {
     fdrfac = cut(x$fdr, c(-.01, .05, .1, .25, .5, 1.01))
     relevel(fdrfac, "(0.5,1.01]")
     },
    mafcats = function(x) {
     maffac = cut(x$MAF,c(-0.01,.05, .1, .25, .51))
     relevel(maffac, "(-0.01,0.05]")
     },
    caddcats = function(x){
     cut(x$PHRED, c(-.01, 5, seq(10, 30, 10 ), 60))
    }
    ),
    cutts = c(-0.01,seq(0.015,.12,.015),.15), 
    names2check=GGtools:::.standardNames, maxit=30, savePinfer=FALSE
   ) {  # finish list and function arg paren

requireNamespace("foreach")

#
# subroutines
#

.discretize_dt = function( dtab, txlist ) {
   newfacs = lapply(txlist, function(f) f(dtab) )
   newvar = names(txlist)
   for (i in 1:length(newvar))
     dtab[[newvar[i]]] = newfacs[[i]]
   dtab
   }  # end .discretize_dt

.redu.fdr = function(dtab) {
 setkey(dtab,snp,score)  # now sorted
 inds = cumsum(dtab[,.N,by="snp"]$N)
 sco = as.numeric(dtab[inds,score])
 pnames = grep("permScore", names(dtab), value=TRUE)
 rdtab = dtab[inds,]
 perms = lapply(pnames, function(x) {
   setkeyv(dtab, c("snp", x))  # resort within SNP by successive perm scores
   inds = cumsum(dtab[,.N,by="snp"]$N)  # take max per snp
   dtab[inds,x,with=FALSE][[x]]  # retain the max
   })
 rdtab$fdr = pifdr(sco, unlist(perms))
 rdtab
}  # end .redu.fdr

.discmods = function( dtab, prefix, folder,
   discfmlas = discfmlas_in, cutts, inmaxit=30) {
  requireNamespace("foreach")
#  curwd = getwd()
#  if (!file.exists(folder)) dir.create(folder)
#  setwd(folder)
#  on.exit(setwd(curwd))
  intrain = 1*( dtab$seqnames %in% paste0("chr", seq(1,22,2)))
  dtab$fdrlt10 = 1*(dtab$fdr < 0.10)
  dtab$fdrlt05 = 1*(dtab$fdr < 0.05)
  dtab$fdrlt01 = 1*(dtab$fdr < 0.01)
  dtab$fdrlt005 = 1*(dtab$fdr < 0.005)
  train = dtab[which(intrain==1),]
  test = dtab[-which(intrain==1),]

  outs = foreach(i=1:length(discfmlas)) %dopar% {
     tmp = bigglm( discfmlas[[i]], data=train, family=binomial(), chunksize=50000, maxit=inmaxit )
     tpreds = predict(tmp, newdata=test, type="response" )
     rocpreds = try(prediction( tpreds, labels=test$isgwashit ))
     AUC = NA
     if (!inherits(rocpreds, "try-error")) AUC = performance(rocpreds, "auc")
     infer = bigglm( discfmlas[[i]], data=dtab, family=binomial(), chunksize=50000, maxit=inmaxit )
     pinfer = NULL
     if (savePinfer) pinfer = predict(infer, newdata=dtab, type="response")
     ans = list(coefs=coef(tmp), auc=AUC, mat=summary(tmp)$mat, infcoefs=coef(infer),
      infvcov=vcov(infer), infmat=summary(infer)$mat, pinfer=pinfer)
    
  }
  names(outs) = names(discfmlas)
  dfobn = paste0(prefix, "_discfmlas")
  assign(dfobn, discfmlas)
 # save( list=dfobn, file=paste0(dfobn, ".rda"))
  oobn = paste0(prefix, "_outs")
  assign(oobn, outs)
 # save( list=oobn, file=paste0(oobn, ".rda"))

ns = names(discfmlas)
pns = paste0("p", ns)

test$fdrlt10 = 1*(test$fdr < .10)

tabs = coeflist = list()
for (i in 1:length(discfmlas)) {
 cat(i)
 mm = model.matrix( discfmlas[[i]], data=test )
 test[[ pns[i] ]] <-  curp <- plogis(mm %*% outs[[i]][[1]])
 coeflist[[ pns[i] ]] = summary(biglm( isgwashit ~ cut(curp,cutts)-1, data=test ))$mat
 tabs[[ pns[i] ]] = table(cut(curp,cutts))
 }
tobn = paste0(prefix, "_test")
assign(tobn, test)
#save(list=tobn, file=paste0(tobn, ".rda"))
cobn = paste0(prefix, "_coeflist")
assign(cobn, coeflist)
#save(list=cobn, file=paste0(cobn, ".rda"))
tabobn = paste0(prefix, "_tabs")
assign(tabobn, tabs)
#save(list=tabobn, file=paste0(tabobn, ".rda"))
list(coeflist=coeflist, tabs=tabs, outs=outs, dimdtab=dim(dtab), dimtest=dim(test), dimtrain=dim(train))
}  # end .discmods


.valdt = function(x) all(names2check %in% names(x))

#
# end "subroutines"
#

#
# execute the appraisal, if wanted
#
   if (!is.null(names2check)) stopifnot( .valdt(dtab) )

   if (discretize) {
    assign(obn1 <- paste0(prefix, "_dt"), .discretize_dt(dtab, txlist))
#    save(list=obn1, file=paste0(obn1, ".rda"))
    disctab = get(obn1)
    }
   else disctab = dtab
   
   if (reduceToSNP) {
    obn2 = paste0(prefix, "bySNP_dt")
    assign(obn2, .redu.fdr(disctab))
#    save(list=obn2, file=paste0(obn2, ".rda"))
    discBySnp = get(obn2)
    }
   else discBySnp = dtab
   
   .discmods( discBySnp, prefix=prefix, folder=folder, cutts=cutts,
     inmaxit=maxit )
}


calfig = function (colist, tabs, ind = 10, 
    hfudgetxt = 0.0155, tickend=.16, tickgap=.02, ylimin=c(-.01, .16),
    xlimin = c(-.01, .16), fraccex=.8, fuselast=0) 
{
#
# subroutine that assumes a 'cut' is used in the
# calibration assessment, it pulls out numeric boundaries
# of cut intervals
#
    getmidpts = function(coefrn) {
        lima = gsub(".*\\(", "", coefrn)
        limb = gsub("]", "", lima)
        limc = strsplit(limb, ",")
        limd = sapply(limc, as.numeric)
        apply(limd, 2, mean)
    }
#
# end subroutine
#
    midcuts = getmidpts(rownames(colist[[ind]]))
    plot(0, 0, ylim = ylimin, xlim = xlimin,
        pch = " ", xlab = "predicted", ylab = "empirical", axes = FALSE)
    axis(1, at = seq(0, tickend, tickgap))
    axis(2, at = seq(0, tickend, tickgap))
    fracs = paste(nums <- round(tabs[[ind]] * colist[[ind]][, 1]), tabs[[ind]], 
        sep = "/")
    if (fuselast == 0) {
    points(midcuts, colist[[ind]][, 1])
    text(midcuts + hfudgetxt, colist[[ind]][, 1], labels = fracs, 
        cex = fraccex)
    }
    else {
         l.m = length(midcuts)
         idrop = (l.m-fuselast+1):l.m
         midcuts = c(midcuts[-idrop], mean(midcuts[idrop]))
         ndrp = tabs[[ind]][idrop]
         cdrp = colist[[ind]][idrop,1]
         hdrp = cdrp*ndrp
         newfrac = sum(hdrp)/sum(ndrp)
         newfracc = paste(round(sum(hdrp),0), sum(ndrp), sep="/")
         fracs = c(fracs[-idrop], newfracc)
         newco = colist[[ind]][-idrop,1]
         newco = c(newco, newfrac)
         points(midcuts, newco)
         text(midcuts+hfudgetxt, newco, labels=fracs, cex=fraccex)
         }
    abline(0, 1, lty = 2, col = "gray")
}

dtToAppr = function(dt, prunevec, prefix, fmlas,
   txlist) {
 pruned = dt[ which(dt$snp %in% prunevec), ]
 appraise( pruned, discretize=TRUE, reduceToSNP=TRUE, prefix=prefix,
   discfmlas_in = fmlas, names2check=NULL, txlist=txlist )
 appraise( dt, discretize=TRUE, reduceToSNP=TRUE, prefix=paste0(prefix, "_UNP"),
   discfmlas_in = fmlas, names2check=NULL, txlist=txlist )
}

sensFromDT = function(dt,
    targfdrs = c(.1, .05, .01), parmslist =
    list(mafs=c(.025, .05, .075), dists=c(10000, 50000, 100000, 200000))) {
  snpenum = eqsens_dt(dt,
    by = "snps", targfdrs=targfdrs, parmslist=parmslist )
  probenum = eqsens_dt(dt,
    by = "probes", targfdrs=targfdrs, parmslist=parmslist )
  list(snpenum=snpenum, probenum=probenum)
}

#bindgwava = function(gwavadt, eqdt) {
# meq = match(eqdt$snp, gwavadt$snp)
# eqdt$gwava_tss = gwavadt$tss[meq]
# eqdt$gwava_unmat = gwavadt$unmatched[meq]
# eqdt$gwava_regi = gwavadt$region[meq]
# eqdt
#} 

# idea here is that caddPack is on board
addcadd = function(dt, binder=bindcadd) {
 allb = foreach(x=1:22) %dopar% { binder(dt, x) }
 do.call(rbind, allb)
}




waldtests = function(ob, modname="basehmmall", type="pruned",
   mlog10p=TRUE) {
 requireNamespace("aod")
 if (type=="pruned") sel = getPruned(ob)
 else sel = getUnpruned(ob)
 modstr = sel@outs[[modname]]
 infco = modstr$infcoefs
# tags = unique(gsub(".*)\\(", "\\1", names(infco))[-1]) # drop intercept
 tags = unique(gsub("(.*)(878.*)|(\\(.*)", "\\1", names(infco)[-1])) # drop int
 termvecs = lapply(tags, function(x) grep(x, names(infco)))
 wt = lapply(termvecs, function(x) wald.test(modstr$infvcov, infco, Terms=x))
 ans = t(sapply(sapply(wt, function(x) x$result), force))
 rownames(ans) = tags
 if (mlog10p) {
    ans[,"P"] = -log10(ans[,"P"])
    if (!all(isf <- is.finite(ans[,"P"]))) {
          ans[ -which(isf), "P" ] = 16
          colnames(ans)[which(colnames(ans)=="P")] = "-log10P"
          }
    }
 ans
}

procnames = function(covec, dropint=TRUE, dropper=function(x) {
   if (any(x=="none")) which(x=="none") else NA },
   premap = c("distcats"="DIST(bp)", 
        "mafcats"="MAF", "chrom878cat"="CHROMHMM", "caddcats"="CADD", 
        "fdrcats"="eQTL FDR", "gwavucat"="GWAVA(u)"),
   tags = "mafcats|distcats|chromcat878|caddcats|fdrcats|gwavucat") {
    if (dropint) covec=covec[-1]
    na = names(covec)
    NC = length(covec)
#
# most of the model coefficients have names with cut()-generated mangling
# we want whatever precedes the (
#
    prena = gsub("\\(.*", "", na)
#
# chromcat878 needs different handling
#
    chkchromcat = grep("chromcat", na)
    if (length(chkchromcat) > 0) 
        prena[chkchromcat] = "chrom878cat"
    nrle = rle(prena)$lengths
    nrle = cumsum(c(1, nrle[-length(nrle)]))
#
    cats = gsub(tags, "", na)
    bcats = rep(NA, length(cats))
    bcats[nrle] = premap[unique(prena)]
    dr = dropper(cats)
    if (!is.na(dr)) { 
       bcats=bcats[-dr]
       cats=cats[-dr]
       }
    cbind(bcats,cats)
}

fplot = function(obj, mname, dropper = 
    function(x) { if (length(dr <- grep("none$", x))>0) dr else NA }, ...) {
   ltext = procnames(co <- obj@outs[[mname]]$infcoef, dropper=dropper)
   mat = obj@outs[[mname]]$infmat[-1,]
   if (!is.na(dr <- dropper(names(co)[-1])))
     mat = mat[-dr,]
   forestplot( ltext, mat[,1], mat[,2], mat[,3], ...)
}
