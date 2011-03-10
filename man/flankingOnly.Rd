\name{flankingOnly}
\alias{flankingOnly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
helper function to create GRanges representing 
intervals symmetrically abutting ranges
}
\description{
helper function to create GRanges representing 
intervals symmetrically abutting ranges
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
flankingOnly(gr, rad = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr}{
\code{\link[GenomicRanges]{GRanges-class}} instance
}
  \item{rad}{
`radius' or (common) size of each range that surrounds the given range
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
If \code{gr} is length n, value is a
GRanges instance of length 2n
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
library(GenomicRanges)
dem = GRanges(seqnames="chr1", IRanges(c(1e6,2e6), c(1.5e6, 3e6)))
flankingOnly(dem,1e5)
flankingOnly(dem+1e5,1e5)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
