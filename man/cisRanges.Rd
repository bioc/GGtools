\name{cisRanges}
\alias{cisRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create GRanges instance for intervals cis to a set of genes
}
\description{
create GRanges instance for intervals cis to a set of genes
}
\usage{
cisRanges(probeids, chr, anno, radius = 5e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probeids}{
character vector of array probe identifiers for annotation
mapped by package identified in \code{anno}
}
  \item{chr}{
a string to be used for \code{seqnames} component of
GRanges result
}
  \item{anno}{
string identifying an annotation package from which
locations will be derived
}
  \item{radius}{
count of basepairs upstream and downstream of (first) CHRLOC
entry defining the range
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a \code{\link[GenomicRanges]{GRanges-class}} instance
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
library(illuminaHumanv1.db)
g20 = get("20", revmap(illuminaHumanv1CHR))
m1 = intersect(g20, mappedkeys(illuminaHumanv1CHRLOC))
cisRanges(m1[1:10], "chr20", "illuminaHumanv1.db")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
