\name{bestCis}
\alias{bestCis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
extract best (or all) cis-associated eQTL from a multffmgr instance
}
\description{
extract best (or all) cis-associated eQTL from a multffmgr instance
}
\usage{
bestCis(ffmgr, slranges, radius = 1e+06, ffind = 1, anno, ncores = 10)
#allCisP_1sided(ffmgr, slranges, radius = 1e+06, ffind = 1, anno, ncores = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ffmgr}{
manager object, output of multffCT or diagffCC
}
  \item{slranges}{
snp locations RangedData instance
}
  \item{radius}{
number of bases up and down stream to declare cis
}
  \item{ffind}{
index into fflist component of manager for eQTL associations cores
}
  \item{anno}{
character atom naming annotation package for resolution of colnames of ffmgr matrix
}
  \item{ncores}{
number of cores for mclapply to use
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
for \code{bestCis}, data frame with genes as rows, rsnum and chisq(df) scores,
with df and gene and SNP locations as columns.

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
example(diagffCC)
data(snpLocs20)
bestCis(ff, snpLocs20, anno="illuminaHumanv1.db")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
