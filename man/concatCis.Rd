\name{concatCis}
\alias{concatCis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
combine a list of cisRun instances to a single instance
}
\description{
combine a list of cisRun instances to a single instance,
with ad hoc metadata combination
}
\usage{
concatCis(crl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{crl}{
list of instances of \code{\link{cisRun-class}}
}
}
\details{
the metadata for the output is a list with elements
call and config as required, derived from first
element of the input; the extras component
holds the metadata elements of the remaining
input list elements
}
\value{
a cisRun instance
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
example(All.cis)
concatCis(f1, f1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
