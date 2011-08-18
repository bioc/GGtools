\name{genewiseFDRtab}
\alias{genewiseFDRtab}
\alias{policyFDRtab}
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
	targp = c(0.95, 0.975, 0.99, 0.995), 
	folderstem = "fdrf", geneApply = mclapply)
policyFDRtab(sms, rhs, universe=featureNames(sms),
  	policyClo=function(mgr) function(x)topFeats(probeId(x),
     	mgr=mgr, ffind=1, n=1), nperm=1, seed=1234, 
	targp=c(.95, .975, .99, .995),
        folderstem="fdrf", geneApply=mclapply) 

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
thresholds
}
  \item{folderstem}{
name of folder created will start with this string; 'p' will be prepended
for permutation folders
}
  \item{geneApply}{
iteration (apply-like) function to be used to iterate over genes
}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
