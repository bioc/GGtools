\name{multffCT}
\alias{multffCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
run a parallelized cis-trans eQTL search
}
\description{
run a parallelized cis-trans eQTL search
}
\usage{
multffCT(listOfSms, gfmla, geneinds = 1:10, harmonizeSNPs = FALSE, targdir = ".", runname = "foo", overwriteFF = TRUE, fillNA = TRUE, ncores = 2, vmode = "single", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listOfSms}{
%%     ~~Describe \code{listOfSms} here~~
}
  \item{gfmla}{
%%     ~~Describe \code{gfmla} here~~
}
  \item{geneinds}{
%%     ~~Describe \code{geneinds} here~~
}
  \item{harmonizeSNPs}{
%%     ~~Describe \code{harmonizeSNPs} here~~
}
  \item{targdir}{
%%     ~~Describe \code{targdir} here~~
}
  \item{runname}{
%%     ~~Describe \code{runname} here~~
}
  \item{overwriteFF}{
%%     ~~Describe \code{overwriteFF} here~~
}
  \item{fillNA}{
%%     ~~Describe \code{fillNA} here~~
}
  \item{ncores}{
%%     ~~Describe \code{ncores} here~~
}
  \item{vmode}{
%%     ~~Describe \code{vmode} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (listOfSms, gfmla, geneinds = 1:10, harmonizeSNPs = FALSE, 
    targdir = ".", runname = "foo", overwriteFF = TRUE, fillNA = TRUE, 
    ncores = 2, vmode = "single", ...) 
{
    theCall = match.call()
    require(ff, quietly = TRUE)
    .checkArgsMF(listOfSms, fgmla, geneinds, targdir, runname)
    listOfSms = reduceGenes(listOfSms, geneinds)
    if (harmonizeSNPs) 
        listOfSms = makeCommonSNPs(listOfSms)
    else if (!isTRUE(checkCommonSNPs(listOfSms))) 
        stop("harmonizeSNPs = FALSE but SNPs not common across listOfSms, run makeCommonSNPs")
    sumScores2ff(listOfSms, gfmla, targdir, runname, theCall, 
        overwriteFF = overwriteFF, fillNA = fillNA, write = TRUE, 
        ncores = ncores, vmode = vmode, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
