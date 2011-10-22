\name{pcChooser}
\alias{pcChooser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
utility to assist in choosing number of PCs to remove
owing to expression heterogeneity
}
\description{
utility to assist in choosing number of PCs to remove
owing to expression heterogeneity -- only cis testing as of jan 2011
}
\usage{
pcChooser(sms, cand = c(1, 10, 15, 20, 25, 30, 40), fmla, radius = c(1e+05), chr = "20", smlc = "20", geneApply = lapply, pvals = c(1e-06, 1e-07, 1e-08, 1e-09), ncore=NULL, 
 ffind=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sms}{
instance of \code{\link[GGBase]{smlSet-class}}
}
  \item{cand}{
number of PCs to be excluded in successive runs
}
  \item{fmla}{
formula to be used by \code{\link{cisProxScores}}
}
  \item{radius}{
number of basepairs up and downstream from gene boundaries to be
checked for eQTL
}
  \item{chr}{
chromosome for current run, for use in space selection for
GRanges-associated SNP addressing
}
  \item{smlc}{
name of chromosome in \code{names(smList(sms))} for this run 
}
  \item{geneApply}{
iterator to be used for genes
}
  \item{pvals}{
upper bounds on p-values to declare eQTL present
}
  \item{ncore}{if set to numeric value, options(cores=ncore) will 
be executed by this function, useful if geneApply=mclapply}
  \item{ffind}{chrom selector passed to cisProxScores, typically default is
appropriate choice}
  \item{\dots}{
passed to \code{\link{cisProxScores}}
}
}
\details{
The idea is that we want to maximize the number of eQTL declared, and
that there will be diminishing returns as the number of PCs included
grows.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
matrix with columns corresponding to \code{cands} and 
rows corresponding to \code{pvals} -- the row names are
the chi-squared threshold values for \code{\link[snpStats]{snp.rhs.tests}}
results
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
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
if (!exists("hmceuB36.2021")) hmceuB36.2021 <- getSS("GGtools", c("20", "21"))
library(illuminaHumanv1.db)
g20 = get("20", revmap(illuminaHumanv1CHR))
g20 = intersect(g20, featureNames(hmceuB36.2021))[1:40]
pcChooser( hmceuB36.2021[probeId(g20),], cand=c(7,9,11), fmla=~male,
  radius=1e6, chr="20", smlc="20", geneApply=lapply, pvals=10^(-c(3:5)))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
