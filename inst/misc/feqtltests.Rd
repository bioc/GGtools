\name{feqtltests}
\alias{feqtltests}
\title{
compute and manage eqtl tests based on out-of-memory representation
}
\description{
compute and manage eqtl tests based on out-of-memory representation
}
\usage{
feqtltests(egtSet, rhs = ~1, nchunk = 10, runname = "foo", targdir = "foo", geneApply = lapply, chromApply = lapply, shortfac = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egtSet}{
instance of \code{\link{egtSet-class}}
}
  \item{rhs}{
formula, with no dependent variable, with variables taken from pData(egtSet)
}
  \item{nchunk}{
number of pieces into which testing for each chromosome is broken
}
  \item{runname}{
token for identifying output files, used in prefix
}
  \item{targdir}{
folder where ff files are to be written
}
  \item{geneApply}{
lapply-like function to be used while iterating over genes
}
  \item{chromApply}{
lapply-like function to be used while iterating over chromosomes
}
  \item{shortfac}{
scaling factor -- we are storing results as short ints, and will upscale by this factor to preserve
come fractional value of each score
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
instance of \code{\link{eqtlTestsManager-class}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
experimental
}

%i ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
example("egtSet-class")  # constructs fceu
library(illuminaHumanv1.db)
cptag = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
indc = which(featureNames(fceu) == cptag[1])
#
# get a set of additional genes on chr20
all20 = get("20", revmap(illuminaHumanv1CHR))
g20 = unique(c(all20[1:10], cptag))
lkf = feqtltests( fceu[g20,], ~male )
lkf
f2 = list(fceu[g20,], fceu[g20,])
m2 = mfeqtltests(f2, list(~1,~1), folder="abc", runname="def")
m2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
