\name{geneTrack}
\alias{geneTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a RangedData structure with multffCT test results (as -log10 p values by default)
}
\description{
create a RangedData structure with multffCT test results (as -log10 p values by default)
}
\usage{
geneTrack(mgr, gn, chrtag, locdata, dropDups = TRUE, mlog10p = TRUE, minchisq = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mgr}{
an instance of \code{\link{multffManager-class}}.
}
  \item{gn}{
a character string naming a `gene' (typically name for a microarray probe)
}
  \item{chrtag}{
the name of the chromosome for which SNP scores are desired, names of \code{mgr[["fflist"]]}
}
  \item{locdata}{
a RangedData instance with ranges defining SNP locations (0 width) and names giving rs numbers
or any other SNP identifiers indexing rows in \code{mgr[["fflist"]]} ff matrices.
}
  \item{dropDups}{
logical: should duplicated SNP regions be dropped?
}
  \item{mlog10p}{
logical: should the score generated be -10 log p (if FALSE, the chi-squared variate with
\code{mgf[["df"]]} degrees of freedom is used)
}
  \item{minchisq}{
ignore
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The structure provided as \code{locdata} is filtered for SNP that are
tested in \code{mgr} and scores are added in \code{score} element

export using rtracklayer to visualize series of scores on genomic coordinates

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# runs interactively but not in check on windows
if (.Platform$OS.type != "windows") {
example(multffCT)
dems
g1 = colnames(dems$fflist[[1]])[1]
data(snpLocs_21)
sco = geneTrack( dems, g1, "21", snpLocs_21 )
sco
library(rtracklayer)
export(sco, con=paste(g1, ".wig", sep=""))
readLines(paste(g1, ".wig", sep=""), n=10)
#
# now add to genome browser as a custom track
#
# if you want to modify aspects of the display as a track, use, e.g.,
# nsco = as(sco, "UCSCData")
# nsco@trackLine@name = "[genename]"  etc.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
