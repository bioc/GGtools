setMethod("snps", "racExSet", function(x) get("racs",x@racAssays))
setMethod("exprs", "racExSet", function(object) get("exprs",object@assayData))
setMethod("racAssays", "racExSet", function(x) x@racAssays)
setMethod("snpNames", "racExSet", function(x) featureNames(x@racAssays))
setMethod("rarebase", "racExSet", function(x) x@rarebase)
setMethod("SNPalleles", "racExSet", function(x) x@SNPalleles)



setMethod("snpScreen", c("racExSet", "snpMeta", "genesym", "formula", "function", "numeric"),
   function (racExSet, snpMeta, gene, formTemplate = ~., fitter = lm, 
      gran) 
  {
      runTemplate = function(x, y) {
          z = as.character(x)
          z[2] = gsub("\\.", y, z[2])
          as.formula(z)
      }
      psid = getpsid(gene, annotation(racExSet))
      y = exprs(racExSet)[psid, ]
      outco = list(y)
      names(outco) = as.character(gene)
      nsnp = length(sn <- snpNames(racExSet))
      snpstodo = sn[inuse <- seq(1, nsnp, gran)]
      if (any(is.na(snpstodo))) snpstodo = snpstodo[!is.na(snpstodo)]
      allpos = get("meta", snpMeta@meta)$pos
      allsn = rownames(get("meta", snpMeta@meta))
      names(allpos) = allsn
      snpstodo = intersect(snpstodo, allsn)
      locs = allpos[snpstodo]
      fittertok = deparse(substitute(fitter))
      callsave = match.call()
      out = list()
      if (fittertok %in% c("fastAGM", "fastHET")) {
        tmp = snps(racExSet)[snpstodo,]
        bad = apply(tmp,1,function(x) any(is.na(x)))
        if (any(bad)) {
           warning("some genotype results had missing values; associated SNPs are dropped completely in this version when fastAGM is used.")
           snpstodo = snpstodo[-which(bad)]
           }
        locs = allpos[snpstodo]
        if (fittertok == "fastAGM") ans = fastAGM(snps(racExSet)[snpstodo,], y)
        else if (fittertok == "fastHET") ans = fastHET(snps(racExSet)[snpstodo,], y)
        return(new("snpScreenResult", call=callsave, locs=locs, 
            chr=chromosome(snpMeta), fittertok=fittertok, gene=as.character(gene), ans))
      }
      for (i in 1:length(snpstodo)) {
          if (options()$verbose == TRUE) {
              if (i%%100 == 0) 
                  cat(i)
          }
          fm = runTemplate(formTemplate, snpstodo[i])
          out[[i]] = try(oneFit(racExSet, outco, fm, fitter))
      }
      names(out) = snpstodo
      new("snpScreenResult", call = callsave, locs = locs, chr = chromosome(snpMeta), 
        fittertok = fittertok, gene=as.character(gene), out)
})

setMethod("snpScreen", c("racExSet", "snpMeta", "genesym", "formula", "function", "missing"),
   function (racExSet, snpMeta, gene, formTemplate = ~., fitter = lm, 
      gran) 
  {
#
# absurd to replicate all this code for missing gran, but
# madness of fitter deparsing compels it.  need to introduce
# GGfitter class and then dispatch when it is present, otherwise
# call the function
#
      gran = 1
      runTemplate = function(x, y) {
          z = as.character(x)
          z[2] = gsub("\\.", y, z[2])
          as.formula(z)
      }
      psid = getpsid(gene, annotation(racExSet))
      y = exprs(racExSet)[psid, ]
      outco = list(y)
      names(outco) = as.character(gene)
      nsnp = length(sn <- snpNames(racExSet))
      snpstodo = sn[inuse <- seq(1, nsnp, gran)]
      if (any(is.na(snpstodo))) snpstodo = snpstodo[!is.na(snpstodo)]
      allpos = get("meta", snpMeta@meta)$pos
      allsn = rownames(get("meta", snpMeta@meta))
      names(allpos) = allsn
      snpstodo = intersect(snpstodo, allsn)
      locs = allpos[snpstodo]
      fittertok = deparse(substitute(fitter))
      callsave = match.call()
      out = list()
      if (fittertok %in% c("fastAGM", "fastHET")) {
        tmp = snps(racExSet)[snpstodo,]
        bad = apply(tmp,1,function(x) any(is.na(x)))
        if (any(bad)) {
           warning("some genotype results had missing values; associated SNPs are dropped completely in this version when fastAGM is used.")
           snpstodo = snpstodo[-which(bad)]
           }
        locs = allpos[snpstodo]
        if (fittertok == "fastAGM") ans = fastAGM(snps(racExSet)[snpstodo,], y)
        else if (fittertok == "fastHET") ans = fastHET(snps(racExSet)[snpstodo,], y)
        return(new("snpScreenResult", call=callsave, locs=locs, 
            chr=chromosome(snpMeta), fittertok=fittertok, gene=as.character(gene), ans))
      }
      for (i in 1:length(snpstodo)) {
          if (options()$verbose == TRUE) {
              if (i%%100 == 0) 
                  cat(i)
          }
          fm = runTemplate(formTemplate, snpstodo[i])
          out[[i]] = try(oneFit(racExSet, outco, fm, fitter))
      }
      names(out) = snpstodo
      new("snpScreenResult", call = callsave, locs = locs, chr = chromosome(snpMeta), 
        fittertok = fittertok, gene=as.character(gene), out)
})

setMethod("show", "racExSet", function(object) {
    cat("racExSet instance (SNP rare allele count + expression)\n")
    cat("rare allele count assayData:\n")
  cat("  Storage mode:", storageMode(racAssays(object)), "\n")
  nms <- selectSome(snpNames(object))
  cat("  featureNames:", paste(nms, collapse=", "))
  if ((len <- length(snpNames(object))) > length(nms))
    cat(" (", len, " total)", sep="")
  cat("\n  Dimensions:\n")
  print(Biobase:::assayDataDims(racAssays(object)))
  cat("\nexpression assayData\n")
  cat("  Storage mode:", storageMode(object), "\n")
  nms <- selectSome(featureNames(object))
  cat("  featureNames:", paste(nms, collapse=", "))
  if ((len <- length(featureNames(object))) > length(nms))
    cat(" (", len, " total)", sep="")
  cat("\n  Dimensions:\n")
  print(dims(object))
  cat("\nphenoData\n")
  show(phenoData(object)) 
  cat("\n")
  show(experimentData(object))
  cat("\nAnnotation ")
  show(annotation(object))
    })

make_racExSet = function(exprs, racs, rarebase, SNPalleles, pd, mi, anno) {
    if (!is(exprs, "matrix")) 
        stop("exprs must be of class matrix")
    if (!is(racs, "matrix")) 
        stop("racs must be of class matrix")
    if (!is(pd, "phenoData") & !is(pd, "AnnotatedDataFrame")) 
        stop("pd must be of class phenoData or AnnotatedDataFrame")
    names(SNPalleles) = rownames(racs)
    names(rarebase) = rownames(racs)
    new("racExSet", exprs=exprs, racs=racs, rarebase=rarebase, 
        SNPalleles = SNPalleles,
        phenoData = pd, experimentData = mi, 
        annotation = anno)
}

plot_EvG = function(reset, gene, snpid, anno="hgfocus") {
 gn = getpsid(gene, anno)
 y = exprs(reset)[gn,]
 x = snps(reset)[snpid,]
 plot(x,y,ylab=paste("log", gene, "expression"), xlab=paste("minor allele count,",
  snpid), pch=20, axes=FALSE)
 axis(2)
 axis(1, at=0:2)
}

setGeneric("racAssays<-", function(object,value)standardGeneric("racAssays<-"))
setReplaceMethod("racAssays", c("racExSet", "AssayData"), function(object, value) {
 object@racAssays = value
 object
})

setMethod("[", "racExSet", function(x, i, j, ..., drop=FALSE) {
 if (is(i, "genesym")) callNextMethod()
 else if (is(i, "snpID")) {
    sel = get("racs", x@racAssays)
    sel = sel[i,,drop=FALSE]
    x@racAssays = assayDataNew("lockedEnvironment", racs=sel)
 }
 x
})
