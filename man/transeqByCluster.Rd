\name{transeqByCluster}
\alias{transeqByCluster}
\alias{transeqByChrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convenience functions for trans-eQTL testing
}
\description{
convenience functions for trans-eQTL testing, one assuming
a parallel-based cluster instance is available, one assuming
a chromosome's SNPs will all be candidates for testing
}
\usage{
transeqByCluster(cl, 
  snpchrs = c("chr21", "chr22"), 
  exchrs = 1:22, baseconf, 
  targname = "transrun_", nperm = 1, inseed = 1234, ...)

transeqByChrom(snpchr = "chr22", 
    exchrs = 1:22, baseconf, targname = "transrun_", 
    nperm = 1, inseed = 1234, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{
cluster instance as defined by the parallel package
makeCluster API
}
  \item{snpchrs}{
character vector of tokens to be used to enumerate chromosomes
harboring SNP for testing
}
  \item{snpchr}{
character atom, for transeqByChrom, the chromosome on which
testing will be conducted 
}
  \item{exchrs}{
enumeration of chromosomes harboring expression measures to
be checked for trans association with SNPs
}
  \item{baseconf}{
an instance of \code{\link{TransConfig-class}}
}
  \item{targname}{
folder where scratch results are computed
}
  \item{nperm}{
number of permutations to be used for plug-in FDR
}
  \item{inseed}{
seed to be set before permutations are attempted,
in conjunction with
    RNGkind("L'Ecuyer-CMRG")
}
  \item{\dots}{
not used
}
}
\details{
the \code{\link{TransConfig-class}} instance determines
most of the details of the testing procedure
}
\value{
a data.frame with test results as chisq, and permScore*
with scores obtained after permuting expression against genotype
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
