\name{diagffCC}
\alias{diagffCC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
perform a 'diagonal' cis eQTL search (only check SNPs chromosomally coresident with genes)
}
\description{
perform a 'diagonal' cis eQTL search (only check SNPs chromosomally coresident with genes)
}
\usage{
diagffCC(sms, gfmla, targdir = ".", runname = "foo", overwriteFF = TRUE, ncores = 2, vmode = "short", shortfac = 100, mc.set.seed = TRUE, fillNA = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sms}{
smlSet
}
  \item{gfmla}{
formula with right-hand side specifying covariates, dependent variable
should be 'gs'
}
  \item{targdir}{
folder to hold results
}
  \item{runname}{
arbitrary distinguishing tag
}
  \item{overwriteFF}{
preserve preexisting FF files if FALSE
}
  \item{ncores}{
number of cores to use with multicore
}
  \item{vmode}{
can be "short" to use efficient space
}
  \item{shortfac}{
amount to scale short ints by to preserve some precision
}
  \item{mc.set.seed}{
as in multicore
}
  \item{fillNA}{
when test cannot be performed (eg due to monomorphy) fill in with chisq(1) if true
}
  \item{\dots}{
passed to snp.rhs.tests of snpMatrix
}
}
\details{
uses annotation package specified in annotation slot of smlSet (which should have
.db suffix) to get list of genes on each chromosome present in smlSet
}
\value{
a multffManager instance
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
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
data(hmceuB36.2021)
library(illuminaHumanv1.db)
g20 = get("20", revmap(illuminaHumanv1CHR))[1:10]
g21 = get("21", revmap(illuminaHumanv1CHR))[1:10]
cpn = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
g20 = c(g20,cpn)
hh = hmceuB36.2021[probeId(c(g20,g21)),]
owd = getwd()
setwd(ind <- tempdir())
print(ind)
ff = diagffCC( hh, gs~male, runname="test")
ff
# we know the following should have a score above 50
ff[ rsid("rs6060535"), probeId(cpn) ]
setwd(owd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
