
setClass("meFiles", representation(exfile="character",
  gtfile="character", outfile="character"))
setGeneric("exfile", function(x)standardGeneric("exfile"))
setGeneric("gtfile", function(x)standardGeneric("gtfile"))
setGeneric("outfile", function(x)standardGeneric("outfile"))
setMethod("exfile", "meFiles", function(x)x@exfile)
setMethod("gtfile", "meFiles", function(x)x@gtfile)
setMethod("outfile", "meFiles", function(x)x@outfile)

textAndSlice = function(sms, cind=1) {
 ex = exprs(sms)
 gt = t(smList(sms)[[cind]])
 exfile = tempfile()
 gtfile = tempfile()
 outfile = tempfile()
 write.table(ex, exfile, sep="\t")
 write.table(as(gt, "numeric"), gtfile, sep="\t")
 new("meFiles", exfile=exfile, gtfile=gtfile, outfile=outfile)
}

getAsSlicedData = function(sms, SNPtarget, cind=1) {
 if (!is(SNPtarget, "SlicedData")) stop("target must be SlicedData instance")
 ex = exprs(sms)
 gt = as(t(smList(sms)[[cind]]), "numeric")
 target$CreateFromMatrix(gt)
 new("meFiles", exfile=exfile, gtfile=gtfile, outfile=outfile)
}


.slicedDataDefaults = list(
  fileDelimiter = "\t",
  fileOmitCharacters = "NA",
  fileSkipRows = 1,
  fileSkipColumns = 1,
  fileSliceSize = 2000)

eqtlTests.meText = function( smlSet, runname="20", targdir="cisScratch.me",
      geneApply=lapply, shortfac = 100, checkValid = TRUE, useUncertain= TRUE,
      glmfamily = "gaussian", scoretx = abs,
      matrixEQTL.engine.control = list(output_file_name = outfile(mefob), 
         pvOutputThreshold=1e-3,
         useModel=modelLINEAR, errorCovariance=numeric(), verbose=FALSE, pvalue.hist=FALSE),
      snpSlicedData.control=.slicedDataDefaults,
      geneSlicedData.control=.slicedDataDefaults,
      covarSlicedData.control=.slicedDataDefaults,
      covariates_file_name = character() ) {
#
# this is a very preliminary interface to MatrixEQTL testing procedure
# for performance and inference comparisons
#
# the premise here was that keeping the big data in text files
# was advantageous for reducing memory consumption but this was
# a misconception.  current engine behavior assumes all slices in memory
# -- the eqtlTests.me below does not use text interface but
# directly converts matrices
#
  require(MatrixEQTL)
  thecall = match.call()
  cat("textualizing...")
  mefob = textAndSlice( smlSet )
  cat("done.\n")

  useModel = modelLINEAR; # modelANOVA or modelLINEAR

  snps = SlicedData$new();
  do.call(snps$initFields, snpSlicedData.control )
  snps$LoadFile(gtfile(mefob))

## Load gene expression data

  gene = SlicedData$new();
  do.call(gene$initFields, geneSlicedData.control )
  gene$LoadFile(exfile(mefob))

## Load covariates

  cvrt = SlicedData$new();
  do.call(cvrt$initFields, covarSlicedData.control )
  if(length(covariates_file_name)>0) {
  	cvrt$LoadFile(covariates_file_name);
  }

## Run the analysis

controls = list(matrix=matrixEQTL.engine.control,
 geneSlice=geneSlicedData.control, snpSlice=snpSlicedData.control,
 covarSlice=covarSlicedData.control)
me = do.call(Matrix_eQTL_engine, c(list(snps, gene, cvrt), 
      matrixEQTL.engine.control) )
 
geneNames = featureNames(smlSet)

snpids = colnames(smList(smlSet)[[1]])

ngenes = length(geneNames)
nsnps = length(snpids)

unlink(targdir, recursive=TRUE)
dir.create(targdir)
targff = paste(targdir, "/", runname, ".ff", sep="")

store = ff( initdata=0,
        dim=c(nsnps, ngenes),
        dimnames=list(snpids, geneNames),
        vmode="short",
        filename = targff )

cat("populating store...")

#Browse[1]> x[1:4,]
#       snps         gene statistic       pvalue          FDR      beta
#1  rs407257 GI_4504184-S -19.20425 3.183801e-33 3.078519e-25 -1.326829
#2 rs5760176 GI_4504184-S -15.20243 2.400666e-26 1.160641e-18 -1.208819
#3 rs4822466 GI_4504184-S -13.57040 2.730705e-23 5.824989e-16 -1.227329
#4 rs4822450 GI_4504184-S  13.47201 4.216938e-23 5.824989e-16  1.215933

rowi = match(me$all$eqtls[,"snps"], rownames(store))
coli = match(me$all$eqtls[, "gene"], colnames(store))
store[ cbind(rowi,coli) ] =
   scoretx ( me$all$eqtls[, "statistic"] * shortfac )
cat("done.\n")

new("eqtlTestsManager", fffile=store, call=thecall, sess=sessionInfo(),
  exdate=date(), shortfac=shortfac, geneanno=smlSet@annotation, df=1 )
}

eqtlTests.me = function( smlSet, rhs=~1, runname="20", targdir="cisScratch.me", pvot=.5,
      geneApply=lapply, shortfac = 100, checkValid = TRUE, useUncertain= TRUE,
      glmfamily = "gaussian", scoretx = abs,
      matrixEQTL.engine.control = list( output_file_name="/dev/null",
         useModel=modelLINEAR, errorCovariance=numeric(), verbose=FALSE, pvalue.hist=FALSE),
      snpSlicedData.control=.slicedDataDefaults,
      geneSlicedData.control=.slicedDataDefaults,
      covarSlicedData.control=.slicedDataDefaults,
      covariates_file_name = character() ) {
#
# this is a very preliminary interface to MatrixEQTL testing procedure
# for performance and inference comparisons
#
  require(MatrixEQTL)
  matrixEQTL.engine.control$pvOutputThreshold = pvot
#
  if (length(smList(smlSet)) != 1) stop("please supply smlSet with smList of length 1")
  require(MatrixEQTL)
  thecall = match.call()
#  cat("converting SNP data...")
  snps = SlicedData$new();
  do.call(snps$initFields, snpSlicedData.control )
  snps$CreateFromMatrix( as(t(smList(smlSet)[[1]]), "numeric" ) )

#  useModel = modelLINEAR; # modelANOVA or modelLINEAR


## Load gene expression data

  gene = SlicedData$new();
#  cat("converting expression data...")
  do.call(gene$initFields, geneSlicedData.control )
  gene$CreateFromMatrix(exprs(smlSet))

## Load covariates

  cvrt = SlicedData$new();
  do.call(cvrt$initFields, covarSlicedData.control )
#  cat("converting covariate data...")
  if(!identical(rhs, ~1)){
        dmat = try(model.matrix(rhs, data=pData(smlSet)))
	if (inherits(dmat, "try-error")) stop("rhs could not define model.matrix from pData(smlSet) without error")
  	cvrt$CreateFromMatrix(t(dmat))
  }

# retain call conditions

controls = list(matrix=matrixEQTL.engine.control,
 geneSlice=geneSlicedData.control, snpSlice=snpSlicedData.control,
 covarSlice=covarSlicedData.control)

## Run the analysis

me = do.call(Matrix_eQTL_engine, c(list(snps, gene, cvrt), 
      matrixEQTL.engine.control) )
 
geneNames = featureNames(smlSet)

snpids = colnames(smList(smlSet)[[1]])

ngenes = length(geneNames)
nsnps = length(snpids)

unlink(targdir, recursive=TRUE)
dir.create(targdir)
targff = paste(targdir, "/", runname, ".ff", sep="")

store = ff( initdata=0,
        dim=c(nsnps, ngenes),
        dimnames=list(snpids, geneNames),
        vmode="short",
        filename = targff )

cat("populating store...")

rowi = match(me$all$eqtls[,"snps"], rownames(store))
coli = match(me$all$eqtls[, "gene"], colnames(store))
store[ cbind(rowi,coli) ] =
   scoretx ( me$all$eqtls[, "statistic"] * shortfac )

cat("done.\n")

new("eqtlTestsManager", fffile=store, call=thecall, sess=sessionInfo(),
  exdate=date(), shortfac=shortfac, geneanno=smlSet@annotation, df=1,
  summaryList = controls  )
}
