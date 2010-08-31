\name{getNamedLocs}
\alias{getNamedLocs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get a named vector of SNP locations
}
\description{
get a named vector of SNP locations
}
\usage{
getNamedLocs(slpack = "SNPlocs.Hsapiens.dbSNP.20100427", chrtok)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{slpack}{
location package name
}
  \item{chrtok}{
string giving the token used to obtain a chromosome's worth of SNP locations
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%%  ~Describe the value returned
%%%  If it is a LIST, use
%%%  \item{comp1 }{Description of 'comp1'}
%%%  \item{comp2 }{Description of 'comp2'}
%%% ...
%}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
should eventually give way to GRanges processing
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
nn = getNamedLocs(chr="ch20")
length(nn)
nn[1:5]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
