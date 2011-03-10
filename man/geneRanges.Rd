\name{geneRanges}
\alias{geneRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
construct a RangedData instance for genes enumerated according to an
annotation .db package
}
\description{
construct a RangedData instance for genes enumerated according to an
annotation .db package
}
\usage{
geneRanges(ids, annopkg, extend = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ids}{
character vector
}
  \item{annopkg}{
package that includes CHR, CHRLOC and CHRLOCEND maps for tokens in \code{ids}
}
  \item{extend}{
atomic number of bases to extend ranges from start upstream and from end downstream
}
}
\details{
if no location is available, start is set to 1 and end is set to 2, regardless
of value of \code{extend}
}
\value{
\code{\link[IRanges]{RangedData-class}} instance
}
%\references{
%%% ~put references to the literature/web site here ~
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
\dontrun{
library(illuminaHumanv1.db)
gg = get(c("CPNE1", "BRCA2"), revmap(illuminaHumanv1SYMBOL))
geneRanges(gg, "illuminaHumanv1.db")
}
}
\keyword{ models }
