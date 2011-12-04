  pcChooser = function(sms, cand=c(1,10,15,20,25,30,40), fmla, radius=c(100000), chr="20", smlc="20", 
    geneApply=lapply, pvals=c(1e-6,1e-7,1e-8,1e-9), ncore=NULL, ffind=1, ...) {
 # strictly Cis testing at this time 7 Jan 2011
   if (!is.null(ncore)) options(cores=ncore,mc.cores=ncore)
   if (length(radius) > 1) stop("radius should have length 1")
   require(sms@annotation, character.only=TRUE)
   anpref = gsub(".db", "", sms@annotation)
   p = get(chr, revmap(get(paste(anpref, "CHR", sep=""))))
   p = intersect(p, featureNames(sms))
   sms = sms[probeId(p),]
   sms = sms[chrnum(smlc),]
   if (prod(dim(smList(sms)[[1]])) == 0) stop("no SNPs for this smlc")
   ans = list()
   for (i in 1:length(cand)) {
    cat(i)
    targdir = tempfile()
    dir.create(targdir)
    if (file.exists(targdir)) unlink(targdir, recursive=TRUE)
    ans[[i]] = cisProxScores( clipPCs(sms, 1:cand[i]), fmla, dradset=radius,
            folder=targdir, runname="pcc", geneApply=geneApply, ffind=ffind, ... )
   }
  thresh = qchisq(1-pvals, 1)
  neq = matrix(NA, nr=length(thresh), nc=length(cand))
  for (i in 1:length(pvals)) 
    for (j in 1:length(cand)) neq[i,j] = sum(unlist(ans[[j]][[1]], recursive=TRUE) > thresh[i])
  colnames(neq) = paste("pc1:", cand, sep="")
  rownames(neq) = as.character(round(thresh,3))
  neq
}
  
