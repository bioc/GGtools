\name{multffCT}
\alias{multffCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
parallelized multipopulation cis-trans eQTL searches
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
list of \code{\link{smlSet-class}} instances
}
  \item{gfmla}{
formula with dummy dependent variable and variables on right-hand side drawn from
pData of listOfSms, to be passed to \code{\link[snpMatrix]{snp.rhs.tests}}
}
  \item{geneinds}{
object inheriting from numeric or \code{\link[GGBase]{probeId-class}} to enumerate
genes for analysis
}
  \item{harmonizeSNPs}{
logical indicating whether to skip the call to \code{\link{makeCommonSNPs}} for the
\code{listOfSms}
}
  \item{targdir}{
path to location where ff files will be written
}
  \item{runname}{
tag to be used in ff filenames and for ultimate control object to be serialized
}
  \item{overwriteFF}{
logical indicating whether preexisting ff files with names to be used in this
run should be overwritten (by default they are)
}
  \item{fillNA}{
logical indicating whether array elements corresponding to missing tests should 
be filled with independent chisquared df 1.  Note that concrete reproducibility of sets
of scores that are randomly generated is not achieved if mc.set.seed=TRUE, which
is the default value.
}
  \item{ncores}{
maximum number of cores to be used by \code{\link[multicore]{mclapply}}
}
  \item{vmode}{
mode for numeric storage in ff files, see \code{\link[ff]{vmode}}.
}
  \item{\dots}{
additional arguments for passage to \code{\link[snpMatrix]{snp.rhs.tests}}
}
}
\details{
function constructs nchrom ff files holding sums of chisquared tests across
smlSets supplied in listOfSms, and serializes metadata about them and the
run in \code{[runname].rda}.
}
\value{
a list for inspection, but key result is side effect of writing ff files
and serializing their metadata
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(smlSet.example)
td = tempdir()
od = getwd()
on.exit(unlink(td))
setwd(td)
dem = multffCT( list(smlSet.example, smlSet.example), gs~male, 1:3 )
dir()
}
\keyword{ models }
