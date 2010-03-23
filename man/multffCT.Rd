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
multffCT(listOfSms, gfmlaList, geneinds = 1:10, harmonizeSNPs = FALSE, targdir = ".", runname = "foo", overwriteFF = TRUE, fillNA = TRUE, 
   ncores = 2, mc.set.seed=TRUE, vmode = "single", shortfac=100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listOfSms}{
list of \code{\link{smlSet-class}} instances
}
  \item{gfmlaList}{
list of formulas (associated one to one with components of
\code{listOfSms}) with dummy dependent variable and variables on 
right-hand side drawn from
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
  \item{mc.set.seed}{
as passed to \code{\link[multicore]{mclapply}}}
  \item{vmode}{
mode for numeric storage in ff files, see \code{\link[ff]{vmode}}.  If you
use "short", the "shortfac" will multiply the chisquares so that integer
storage retains some precision (if shortfac = 100, you have two 
digits beyond the decimal point; the short can only represent
0-32767.)  More infrastructure is needed for downstream handling
of the short representation, but it seems worthwhile.
}
  \item{shortfac}{quantity by which short ints will be inflated for storage to allow
more precision in usage}
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
# runs interactively but not in check on windows
if (.Platform$OS.type != "windows") {
 data(smlSet.example)
 td = tempdir()
 od = getwd()
 on.exit(unlink(td))
 setwd(td)
 set.seed(1234)
 dem = multffCT( list(smlSet.example, smlSet.example), list(gs~male, gs~male), 1:3, runname="dem1" )
 set.seed(1234)
 dems = multffCT( list(smlSet.example, smlSet.example), list(gs~male, gs~male), 
    1:3, vmode="short", shortfac=100, runname="dem2" )
 #
 # note that chisq fillin of missing snps make strict numerical reproducibility
 # nontrivial
 dem
 dems
 dir()
}
}
\keyword{ models }
