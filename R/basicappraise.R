
.discfmlas.demo = list(
  base=isgwashit~distcats+mafcats,
  basehmm=isgwashit~distcats+mafcats+chromcat878,
  allfac=isgwashit~distcats+mafcats+chromcat878+caddcats+fdrcats
)

appraise = function(dtab, discretize=TRUE, 
   reduceToSNP=TRUE, prefix, folder=paste0(prefix, "_APPROUT"), 
   discfmlas_in=GGtools:::.discfmlas.demo, txlist = list(
    distcats = function(x) {
     cut(x$mindist, c(-1, seq(0, 200001, 10000)))
     },
    fdrcats = function(x) {
     fdrfac = cut(x$fdr, c(seq(-.01,.25, .05), .3, .4, .5, 1.01))
     relevel(fdrfac, "(0.5,1.01]")
     },
    mafcats = function(x) {
     maffac = cut(x$MAF,c(-0.01,.05, .1, .15, .2, .25, .3, .35, .4,  .51))
     relevel(maffac, "(-0.01,0.05]")
     },
    caddcats = function(x){
     cut(x$PHRED, c(-.01, 0, 1, 2, 4, 6, 8, 10,seq(20, 35, 5), 60))
    }
    ),
    cutts = c(-0.01,seq(0.015,.12,.015),.15)
   ) {  # finish list and function arg paren

require(foreach)

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
   discfmlas = discfmlas_in, cutts) {
  require(foreach)
  curwd = getwd()
  if (!file.exists(folder)) dir.create(folder)
  setwd(folder)
  on.exit(setwd(curwd))
  intrain = 1*( dtab$seqnames %in% paste0("chr", seq(1,22,2)))
  train = dtab[which(intrain==1),]
  train$fdrlt10 = 1*(train$fdr < 0.10)
  train$fdrlt05 = 1*(train$fdr < 0.05)
  train$fdrlt01 = 1*(train$fdr < 0.01)
  train$fdrlt005 = 1*(train$fdr < 0.005)
  test = dtab[-which(intrain==1),]
  test$fdrlt10 = 1*(test$fdr < 0.10)
  test$fdrlt05 = 1*(test$fdr < 0.05)
  test$fdrlt01 = 1*(test$fdr < 0.01)
  test$fdrlt005 = 1*(test$fdr < 0.005)

  outs = foreach(i=1:length(discfmlas)) %dopar% {
     tmp = bigglm( discfmlas[[i]], data=train, fam=binomial(), chunksize=50000, maxit=15 )
     tpreds = predict(tmp, newdata=test, type="response" )
     rocpreds = try(prediction( tpreds, labels=test$isgwashit ))
     AUC = NA
     if (!inherits(rocpreds, "try-error")) AUC = performance(rocpreds, "auc")
     list(coefs=coef(tmp), auc=AUC)
  }
  names(outs) = names(discfmlas)
  dfobn = paste0(prefix, "_discfmlas")
  assign(dfobn, discfmlas)
  save( list=dfobn, file=paste0(dfobn, ".rda"))
  oobn = paste0(prefix, "_outs")
  assign(oobn, outs)
  save( list=oobn, file=paste0(oobn, ".rda"))

ns = names(discfmlas)
pns = paste0("p", ns)

mbase = model.matrix( discfmlas$base, data=test )
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
save(list=tobn, file=paste0(tobn, ".rda"))
cobn = paste0(prefix, "_coeflist")
assign(cobn, coeflist)
save(list=cobn, file=paste0(cobn, ".rda"))
tabobn = paste0(prefix, "_tabs")
assign(tabobn, tabs)
save(list=tabobn, file=paste0(tabobn, ".rda"))
NULL
}  # end .discmods

.standardNames = c("seqnames", "start", "end", "width", "strand", "snp", 
 "snplocs", 
"score", "ests", "se", "fdr", "probeid", "MAF", "dist.mid", "mindist", 
"genestart", "geneend", "isgwashit", "chromcat878", "permScore_1", 
"permScore_2", "permScore_3", "PHRED")

.valdt = function(x) all(.standardNames %in% names(x))

#
# execute the appraisal
#
   stopifnot( .valdt(dtab) )

   if (discretize) {
    assign(obn1 <- paste0(prefix, "_dt"), .discretize_dt(dtab, txlist))
    save(list=obn1, file=paste0(obn1, ".rda"))
    disctab = get(obn1)
    }
   else disctab = dtab
   
   if (reduceToSNP) {
    obn2 = paste0(prefix, "bySNP_dt")
    assign(obn2, .redu.fdr(disctab))
    save(list=obn2, file=paste0(obn2, ".rda"))
    discBySnp = get(obn2)
    }
   else discBySnp = dtab
   
   .discmods( discBySnp, prefix=prefix, folder=folder )
}

calfig = function( colist = MOREAPPR_coeflist, tabs = MOREAPPR_tabs,
   ind = 10, hfudgetxt=.0155 ) {
 midcuts = c(0.0075, 0.0225, 0.0375, 0.0525,
    0.0675, 0.0825, 0.0975, 0.1125, 0.135)
 plot(0,0, ylim=c(-.01,.16), xlim=c(-.01,.16),
    pch=" ", xlab="predicted", ylab="empirical", axes=FALSE)
 axis(1, at=seq(0,.16,.02))
 axis(2, at=seq(0,.16,.02))
 points(midcuts, colist[[ind]][,1])
 fracs = paste(round(tabs[[ind]]*colist[[ind]][,1]), tabs[[ind]],
  sep="/")
 text(midcuts+hfudgetxt, colist[[ind]][,1], labels=fracs, cex=.8)
 abline(0,1,lty=2,col="gray")
}

