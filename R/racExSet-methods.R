
setMethod("snps", c("racExSet","missing"), function(x,chr) get("racs",x@racAssays))
setMethod("exprs", "racExSet", function(object) get("exprs",object@assayData))
setMethod("racAssays", "racExSet", function(x) x@racAssays)
setMethod("snpNames", "racExSet", function(x) featureNames(x@racAssays))
setMethod("rarebase", "racExSet", function(x) x@rarebase)
setMethod("SNPalleles", "racExSet", function(x) x@SNPalleles)

fastAGMfitter = new("GGfitter", name="fastAGM", func=fastAGM)
fastHETfitter = new("GGfitter", name="fastHET", func=fastHET)



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
      callsave = match.call()
      out = list()
      if (is(fitter, "GGfitter")) {
        tmp = snps(racExSet)[snpstodo,]
        bad = apply(tmp,1,function(x) any(is.na(x)))
        if (any(bad)) {
           warning("some genotype results had missing values; associated SNPs are dropped completely in this version when fastAGM is used.")
           snpstodo = snpstodo[-which(bad)]
           }
        locs = allpos[snpstodo]
        if (fitter@name == "fastAGM") ans = fastAGM(snps(racExSet)[snpstodo,], y)
        else if (fitter@name == "fastHET") ans = fastHET(snps(racExSet)[snpstodo,], y)
        return(new("snpScreenResult", call=callsave, locs=locs, annotation=annotation(racExSet),
            chr=chromosome(snpMeta), fitter=fitter, gene=as.character(gene), ans))
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
        fitter = fitter, gene=as.character(gene), annotation=
         annotation(racExSet), out)
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
      callsave = match.call()
      out = list()
      if (is(fitter, "GGfitter")) {
        tmp = snps(racExSet)[snpstodo,]
        bad = apply(tmp,1,function(x) any(is.na(x)))
        if (any(bad)) {
           warning("some genotype results had missing values; associated SNPs are dropped completely in this version when fastAGM is used.")
           snpstodo = snpstodo[-which(bad)]
           }
        locs = allpos[snpstodo]
        if (fitter@name == "fastAGM") ans = fastAGM(snps(racExSet)[snpstodo,], y)
        else if (fitter@name == "fastHET") ans = fastHET(snps(racExSet)[snpstodo,], y)
        return(new("snpScreenResult", call=callsave, locs=locs, annotation=
		annotation(racExSet),
            chr=chromosome(snpMeta), fitter=fitter, gene=as.character(gene), ans))
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
      new("snpScreenResult", call = callsave, locs = locs, 
          chr = chromosome(snpMeta),  annotation=annotation(racExSet),
        fitter= fitter, gene=as.character(gene), out)
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

#plot_EvG = function(reset, gene, snpid, anno="hgfocus") {
# gn = getpsid(gene, anno)
# y = exprs(reset)[gn,]
# x = snps(reset)[snpid,]
# plot(jitter(x),y,ylab=paste("log", gene, "expression"), xlab=paste("minor allele count,",
#  snpid), pch=20, axes=FALSE, cex=1.5)
# axis(2)
# axis(1, at=0:2)
#}

#plot_EvG = function (reset, gene, snpid, anno = "hgfocus",
#     jitfac = 0.5, ...) {
#    gn = getpsid(gene, anno)
#    y = exprs(reset)[gn, ]
#    x = snps(reset)[snpid, ]
#    plot(jitter(x,jitfac), y, ylab = paste("log", gene, "expression"), 
#        xlab = paste("minor allele count,", snpid), pch = 20, 
#        axes = FALSE, cex = 1.9, cex.lab=1.5, ...)
#    axis(2, cex=1.9, cex.axis=1.5)
#    axis(1, at = 0:2, cex=1.9, cex.axis=1.5)
#}

plot_EvG_old = function (reset, gene, snpid, anno = "hgfocus", jitfac = 0.5, 
    ...) 
{
    gn = getpsid(gene, anno)
    y = exprs(reset)[gn, ]
    ry = range(y)
    x = snps(reset)[snpid, ]
    boxplot(split(y, x), axes = FALSE, ylab = paste("log", gene, 
        "expression"), ylim = c(ry[1], ry[2]), xlab = paste("minor allele count,", 
        snpid), range = 0, border = "darkgray", cex.lab = 1.5, names=sort(unique(x)))
    yt = round(ry, 1)
    aty = seq(yt[1], yt[2], 0.2)
    axis(2, cex = 1.9, cex.axis = 1.1, at = aty, las = 2)
    points(jitter(x + 1, jitfac), y, pch = 19)
    axis(1, at = 1:length(unique(x)), cex = 1.9, cex.axis = 1.5, labels=sort(unique(x)))
}



setGeneric("racAssays<-", function(object,value)standardGeneric("racAssays<-"))
setReplaceMethod("racAssays", c("racExSet", "AssayData"), function(object, value) {
 object@racAssays = value
 object
})

setMethod("[", "racExSet", function(x, i, j, ..., drop=FALSE) {
 if (is(i, "genesym")) { #callNextMethod()
    ind = getpsid(i, annotation(x))
    x@assayData = assayDataNew("lockedEnvironment", exprs=exprs(x)[ind,,drop=FALSE])
    }
 else if (is(i, "snpID")) {
    sel = get("racs", x@racAssays)
    sel = sel[i,,drop=FALSE]
    x@racAssays = assayDataNew("lockedEnvironment", racs=sel)
    }
 else if (is(i, "exFeatID")) {
    ind = which( featureNames(x) %in% i )
    x@assayData = assayDataNew("lockedEnvironment", exprs=exprs(x)[ind,,drop=FALSE])
    }
 else stop("selection index must be of class genesym or snpID or exFeatID")
 return(x)
})

setMethod("snpScreen", c("racExSet", "snpMeta", "genesym", "formula", "GGfitter"),
   function (racExSet, snpMeta, gene, formTemplate = ~., fitter = fastAGMfitter, 
      gran) 
  {
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
      callsave = match.call()
      out = list()
        tmp = snps(racExSet)[snpstodo,]
        bad = apply(tmp,1,function(x) any(is.na(x)))
        if (any(bad)) {
           warning("some genotype results had missing values; associated SNPs are dropped completely in this version when fastAGM is used.")
           snpstodo = snpstodo[-which(bad)]
           }
        locs = allpos[snpstodo]
        ans = fitter@func(snps(racExSet)[snpstodo,], y)
        return(new("snpScreenResult", call=callsave, locs=locs, fitter=fitter,
            chr=chromosome(snpMeta), gene=as.character(gene),
		annotation=annotation(racExSet), ans))
      })

setMethod("twSnpScreen", c("racExSet", "snpMeta", "formula", "GGfitter"),
   function (racExSet, snpMeta, formTemplate = ~., fitter = fastAGMfitter){
# transcriptome-wide SnpScreen
  fn = featureNames(racExSet)
  syms = na.omit(unlist(lookUp(fn, annotation(racExSet), "SYMBOL")))
  if (any(dd <- duplicated(syms))) syms = syms[-which(dd)]
  out = list()
  for (i in 1:length(syms))
      out[[i]] = try(snpScreen(racExSet, snpMeta, genesym(syms[i]),
         formTemplate, fitter))
  names(out) = syms
  new("twSnpScreenResult", call=match.call(),
     genes=syms, locs=pos(snpMeta)[rownames(snps(racExSet))], 
        fitter=fitter, out)
  })

setGeneric("imputeSNPFixed", function(x,const=0) standardGeneric("imputeSNPFixed"))
setMethod("imputeSNPFixed", c("racExSet", "missing"), function(x,const=0) {
 ss = snps(x)
 if (!any(is.na(ss))) stop("all SNP non-missing; no imputation")
 ss[is.na(ss)] = const
 x@racAssays = assayDataNew("lockedEnvironment", racs=ss)
 x
})
