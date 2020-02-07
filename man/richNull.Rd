\name{richNull}
\alias{richNull}
\alias{meta.richNull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bind metadata concerning SNP allele frequency and
other aspects of optimized cis-eQTL association to an
mcwBestCis instance
}
\description{
bind metadata concerning SNP allele frequency and
other aspects of optimized cis-eQTL association to an
mcwBestCis instance, to allow conditional FDR computation

}
\usage{
richNull(..., MAFlb = 0.01, npc = 10, radius = 250000, nperm = 1, 
   innerFilt = function(x) x, outerFilt = function(x) x)

meta.richNull(..., MAFlb=.01, npc=10, radius=250000,
   nperm=1, innerFilt=function(x)x, outerFilt=function(x)x) 
#
# internally:
#
#  bigfilt = function(z) 
#    outerFilt(MAFfilter(clipPCs(permEx(innerFilt(z)), 1:npc), lower=MAFlb))
#
}
\arguments{
  \item{\dots}{
should provide bindings for smpack and chrnames,
which will be used to obtain gene/probe locations; see
\code{\link[GGBase:make_smlSet]{getSS}} for information on \code{smpack}
settings.

meta.richNull allows a vector of smpack values bound to \code{smpackvec}
}
  \item{MAFlb}{
lower bound on SNP MAF for null distribution evaluation
}
  \item{npc}{
number of expression principal components to be removed
}
  \item{radius}{
radius used for testing
}
  \item{nperm}{
This establishes how many permutations of expression against genotype will be performed
for this process.
}
  \item{innerFilt}{
function immediately applied to generated smlSet instances 
}
  \item{outerFilt}{
function applied to generated smlSet instances  after clipPCs and MAFfilter are
applied in that order
}
}
\details{
The purpose of \code{richNull} is to obtain realizations from the permutation
distribution of cis-eQTL association statistics, binding information on the
characteristics of the optimal results with the scores.  This allows us to
use conditioning with the realizations from the permutation distribution.
}
\value{
richNull returns a list of \code{nperm} mcwBestCis instances with
additional metadata bound in
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Vince Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
