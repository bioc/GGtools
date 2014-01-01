\name{simpleTiling}
\alias{simpleTiling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a GRanges with a tiling of the human genome
}
\description{
create a GRanges with a tiling of the human genome
}
\usage{
simpleTiling(ntile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ntile}{
%%     ~~Describe \code{ntile} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ntile) 
{
    require(Homo.sapiens)
    hsi = seqinfo(Homo.sapiens)[paste0("chr", 1:22), ]
    GenomicRanges::unlist(tileGenome(hsi, ntile = 100))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
