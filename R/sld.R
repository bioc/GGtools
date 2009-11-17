setGeneric("slimdown", function(x,keepnames) standardGeneric("slimdown"))
setMethod("slimdown", c("gwSnpScreenResult", "logical"), 
    function(x,keepnames=FALSE) {
        cs = lapply(x, "slot", "chisq")
        if (keepnames) {
              ns = lapply(x, "slot", "test.names")
              ns = lapply(ns, function(w) as.integer(gsub("rs", "", w)))
              }
        else ns = lapply(1:length(cs), function(x)NA)
    list(snpnames = ns, chisq=cs)
})

setMethod("slimdown", "multiGwSnpScreenResult", function(x, keepnames=NA) {
     s1 = slimdown(x[[1]], TRUE) # s1$ns is big list of integerized rs numbers
                                 # s1$chisq is list of chisq results for gene 1
     srem = lapply(x[2:length(x)], function(y) slimdown(y, keepnames=FALSE))
     # srem is gene-oriented list with elements having [["chisq"]] 
     # defined for each chrom
     allchi = list()
     for (i in 1:length(s1$chisq)) {
       allchi[[i]] = s1$chisq[[i]]
       for (j in 1:length(srem))
         allchi[[i]] = cbind(allchi[[i]], srem[[j]]$chisq[[i]])
       colnames(allchi[[i]]) = names(x)
       }
     list(snpnames = s1$snpnames, chisq=allchi)
})


# fmla must have gs as dependent variable
slimCisTrans = function( smlSet, genes2do=1:45, ncores=15,
     targdir = "/mnt/data/stvjc/GWAS", fmla=gs~male ) {
  coreinds = rep(1:ncores, each=floor(length(genes2do)/ncores))
  leftover = length(genes2do)-length(coreinds)
  inds = split(genes2do[1:length(coreinds)], coreinds)
  if (leftover > 0) inds[[1]] = c(inds[[1]], genes2do[-(1:length(coreinds))])
  msetup = function(ginds, smlSet, fmla, targdir) {
    require(GSEABase, quietly=TRUE)
    require(annotation(smlSet), character.only=TRUE, quietly=TRUE)
    require(GGtools, quietly=TRUE)
    gs <<- GeneSet(featureNames(smlSet)[ginds],
      geneIdType=AnnotationIdentifier(annotation(smlSet)),
      organism = "Homo sapiens")  # GENERALIZE!
    n = paste(geneIds(gs)[1], "+", sep="")
    assign(n, slimdown(xx <- gwSnpTests(fmla, smlSet)))
    print(xx)
    nf = paste(targdir, paste(n, ".rda", sep=""), sep="/")
    save(list=n, file=nf)
    NULL
   }
  mclapply( inds, msetup, smlSet=smlSet, fmla=fmla, targdir=targdir )
}
    
