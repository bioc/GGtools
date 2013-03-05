\name{All.cis}
\alias{All.cis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function that computes score tests for all SNP cis to genes, with
flexible filtering
}
\description{
function that computes score tests for all SNP cis to genes, with
flexible filtering
}
\usage{
All.cis(smpack, rhs = ~1, nperm = 2, 
  folderstem = "cisScratch", radius = 50000, 
  shortfac = 100, chrnames = "22", smchrpref = "", 
  gchrpref = "", schrpref = "ch", 
  geneApply = lapply, geneannopk = "illuminaHumanv1.db", 
  snpannopk = snplocsDefault(), 
  smFilter = function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.9),  exFilter = function(x) x, keepMapCache = FALSE, 
  SSgen = GGBase::getSS, excludeRadius = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smpack}{
package name for externalized smlSet instance
}
  \item{rhs}{
direct covariate formula for \code{\link[snpMatrix]{snp.rhs.stats}} -- when
using permutation-based FDR it is preferable to work with residuals
}
  \item{nperm}{
number of permutations of expression against genotype for plug-in FDR
computation
}
  \item{folderstem}{
tag to identify a folder for temporary computations
}
  \item{radius}{
number of bases up and downstream of gene location to search for SNP
}
  \item{shortfac}{
factor for scaling short integers to represent association score
}
  \item{chrnames}{
chromosome names to be analyzed
}
  \item{smchrpref}{
prefix to be applied to chrnames elements to select from smpack
}
  \item{gchrpref}{
prefix to be applied to chrnames elements to select gene addresses
}
  \item{schrpref}{
prefix to be applied to chrnames elements for SNPlocs query resolution
}
  \item{geneApply}{
function (like lapply) for iterating over genes, typically will use
mclapply
}
  \item{geneannopk}{
name of package to be used to resolve probe names to gene annotation
}
  \item{snpannopk}{
name of package to be used to resolve SNP identifiers to addresses
}
  \item{smFilter}{
a function that will operate on smlSet instances before testing
}
  \item{exFilter}{
a function that operates on ExpressionSet component of smlSet early on
%%     ~~Describe \code{exFilter} here~~
}
  \item{keepMapCache}{
facility for speeding up the mapping of cis SNP
}
  \item{SSgen}{
special function that can be used to create an smlSet from a nonstandard
package
}
  \item{excludeRadius}{
for binning test procedure
}
  \item{\dots}{
passed to eqtlTests
}
}
\details{
returns score statistics for assocations of all SNP cis to genes,
in a GRanges instance, with range names given by probes; metadata supplied
SNP location, name, and score
}
\value{
GRanges instance
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
  f1 = All.cis("GGdata", chrnames=c("21", "22"))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
