
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


eqtlTests.me = function( smlSet, runname="20", targdir="cisScratch.me",
      geneApply=lapply, shortfac = 100, checkValid = TRUE, useUncertain= TRUE,
      glmfamily = "gaussian", scoretx = abs ) {
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
  snps$fileDelimiter = '\t'; # the TAB character
  snps$fileOmitCharacters = 'NA'; # denote missing values;
  snps$fileSkipRows = 1; # one row of column labels
  snps$fileSkipColumns = 1; # one column of row labels
  snps$fileSliceSize = 2000; # read file in pieces of 2,000 rows
  snps$LoadFile(gtfile(mefob))

## Load gene expression data

  gene = SlicedData$new();
  gene$fileDelimiter = '\t'; # the TAB character
  gene$fileOmitCharacters = 'NA'; # denote missing values;
  gene$fileSkipRows = 1; # one row of column labels
  gene$fileSkipColumns = 1; # one column of row labels
  gene$fileSliceSize = 2000; # read file in pieces of 2,000 rows
  gene$LoadFile(exfile(mefob))

## Load covariates

  covariates_file_name = character()
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = '\t'; # the TAB character
  cvrt$fileOmitCharacters = 'NA'; # denote missing values;
  cvrt$fileSkipRows = 1; # one row of column labels
  cvrt$fileSkipColumns = 1; # one column of row labels
  cvrt$fileSliceSize = 2000; # read file in one piece
  if(length(covariates_file_name)>0) {
  	cvrt$LoadFile(covariates_file_name);
  }

## Run the analysis

me = Matrix_eQTL_engine(
	snps,
	gene,
	cvrt,
	output_file_name = outfile(mefob),
	pvOutputThreshold = .2,
	useModel = modelLINEAR, 
	errorCovariance = numeric(), 
	verbose = FALSE, #TRUE,
	pvalue.hist = 10);

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
