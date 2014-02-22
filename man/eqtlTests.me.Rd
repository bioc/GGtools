\name{eqtlTests.me}
\alias{eqtlTests.me}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
use MatrixEQTL computations and statistics as a back end
to GGtools eqtlTests
}
\description{
use MatrixEQTL computations and statistics as a back end
to GGtools eqtlTests
}
\usage{
eqtlTests.me(smlSet, rhs = ~1, runname = "20", 
  targdir = "cisScratch.me", geneApply = lapply, 
  shortfac = 100, checkValid = TRUE, useUncertain = TRUE, 
  glmfamily = "gaussian", scoretx = abs, 
  matrixEQTL.engine.control = 
    list(output_file_name = tempfile(), 
    pvOutputThreshold = 1e-05, useModel = modelLINEAR, 
    errorCovariance = numeric(), 
    verbose = FALSE, 
    pvalue.hist = FALSE), 
  snpSlicedData.control = .slicedDataDefaults, 
  geneSlicedData.control = .slicedDataDefaults, 
  covarSlicedData.control = .slicedDataDefaults, 
  covariates_file_name = character())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
instance of \code{\link[GGBase]{smlSet-class}}
}
  \item{rhs}{
formula for adjustment of tests for covariates or stratification,
see \code{\link[snpStats]{snp.rhs.tests}}
}
  \item{runname}{
tag used to distinguish emitted files
}
  \item{targdir}{
folder where ff archives will reside
}
  \item{geneApply}{
lapply-like function for iteration over genes,
\code{\link[parallel]{mclapply}} is suitable when in
multicore environments
}
  \item{shortfac}{
scaling factor to increase precision when test results are
stored as short ints in ff
}
  \item{checkValid}{
logical to check validity of input smlSet
}
  \item{useUncertain}{
logical informing snp.rhs.tests that imputed real-valued
B allele counts may be present among genotype data
}
  \item{glmfamily}{
family specification for \code{\link[snpStats]{snp.rhs.tests}}
}
  \item{scoretx}{
function to be applied to MatrixEQTL statistics. Defaults
to abs, for signless association testing 
}
  \item{matrixEQTL.engine.control}{
list of parameters passed to \code{\link[MatrixEQTL]{Matrix_eQTL_engine}}
}
  \item{snpSlicedData.control}{
list of parameters used to define \code{\link[MatrixEQTL]{SlicedData-class}} instances
}
  \item{geneSlicedData.control}{
list of parameters used to define \code{\link[MatrixEQTL]{SlicedData-class}} instances
}
  \item{covarSlicedData.control}{
list of parameters used to define \code{\link[MatrixEQTL]{SlicedData-class}} instances
}
  \item{covariates_file_name}{
if covariates are to be used with MatrixEQTL testing engine, they
reside in this file.  \code{\link[GGBase]{regressOut}} can be used
to avoid this if plug-in FDR are to be used
}
}
\details{
provisional interface
}
\value{
see \code{\link{eqtlTests}}
}
\references{
Shabalin et al Bioinformatics (OUP) 2012
}
%\author{
%%  ~~who you are~~
%}
\note{
intended for simple comparisons
}


%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
g22 = nsFilter( chrFilter( getSS("GGdata", "22"), "22" ), var.cutoff = .8 )
m22 = eqtlTests.me(g22)
}
\keyword{ models }
