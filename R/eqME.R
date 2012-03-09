
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


.slicedDataDefaults = list(
  fileDelimiter = "\t",
  fileOmitCharacters = "NA",
  fileSkipRows = 1,
  fileSkipColumns = 1,
  fileSliceSize = 2000)

eqtlTests.me = function( smlSet, runname="20", targdir="cisScratch.me",
      geneApply=lapply, shortfac = 100, checkValid = TRUE, useUncertain= TRUE,
      glmfamily = "gaussian", scoretx = abs,
      matrixEQTL.engine.control = list(output_file_name = outfile(mefob), pvOutputThreshold=1e-5,
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
store[ cbind(me$all$eqtls[,"snps.id"],  me$all$eqtls[, "gene.id" ]) ] =
   scoretx ( me$all$eqtls[, "statistic"] * shortfac )
cat("done.\n")

new("eqtlTestsManager", fffile=store, call=thecall, sess=sessionInfo(),
  exdate=date(), shortfac=shortfac, geneanno=smlSet@annotation, df=1 )
}
