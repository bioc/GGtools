\name{genewiseFDRtab}
\alias{genewiseFDRtab}
\alias{policywiseFDRtab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
encapsulates testing, permutation and thresholding using permutation
distribution to obtain plug-in estimates of FDR at various
thresholds
}
\description{
encapsulates testing, permutation and thresholding using permutation
distribution to obtain plug-in estimates of FDR at various
thresholds
}
\usage{
genewiseFDRtab(sms, rhs, nperm = 1, seed = 1234, 
	targp = c(0.9, 0.95, 0.975, 0.99, 0.995), 
	folderstem = "fdrf", chromApply = lapply,
	geneApply = lapply, gene2snpList=NULL, ...)
policywiseFDRtab(sms, rhs, nperm = 1, seed = 1234, 
	targp = c(0.95, 0.975, 0.99, 0.995), 
	folderstem = "fdrf", geneApply = lapply, 
  	policyClo=function(mgr) function(x)topFeats(probeId(x),
     	mgr=mgr, ffind=1, n=1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sms}{
instance of \code{\link{smlSet-class}}
}
  \item{rhs}{
formula fragment from tilde; names resolved in \code{pData(sms)}
}
  \item{nperm}{
number of permutations to evaluate
}
  \item{seed}{
\code{set.seed} will be run prior to first permutation
}
  \item{targp}{
percentiles/100 of permutation distribution to be used as
thresholds; not fully integrated into report extraction code as of 10/10/2011,
so it is not advised to manipulate this
}
  \item{folderstem}{
name of folder created will start with this string; 'p' will be prepended
for permutation folders
}
  \item{geneApply}{
iteration (apply-like) function to be used to iterate over genes
}
  \item{gene2snpList}{
list with genes as organizing element, specifying for each
gene which SNP will be retained for testing summaries.  can be
created by the \code{\link{proximityList}} function
}
%  \item{universe}{tests will focus on organizing by the tokens in this vector}
  \item{policyClo}{function accepting an eqtlTests manager, returning function of argument x
that will implement a policy of filtering scores for use in identifying genetic
elements that are associated with expression}
  \item{\dots}{passed to \code{\link{eqtlTests}} as run by genewiseScores}
}
\details{
\code{policyFDRtab} should achieve greater flexibility at the call, allowing
selection policies to be defined over SNP.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(smlSet.example)
applier = lapply
if ("multicore" \%in\% installed.packages()[,1] & .Platform$OS.type != "windows") {
 library(multicore)
 applier = mclapply
}
td = tempdir()
curd = getwd()
setwd(td)
#  this is a poor example but will flag any major faults
f1 = genewiseFDRtab( smlSet.example[1:30,] , ~male, geneApply=applier )
f1
setwd(curd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
