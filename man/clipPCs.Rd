\name{clipPCs}
\alias{clipPCs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simple approach to removal of principal components from smlSet
}
\description{
simple approach to removal of principal components from smlSet
}
\usage{
clipPCs(smlSet, inds2drop, center=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
instance of \code{\link[GGBase]{smlSet-class}}
}
  \item{inds2drop}{
numeric vector of PCs to be eliminated
}
  \item{center}{logical passed to \code{\link[stats]{prcomp}}.}
}
\details{
uses SVD and zeroes out selected eigenvalues before reassembly
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an smlSet instance with transformed expression data
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
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
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(hmceuB36.2021)
hmc = clipPCs(hmceuB36.2021, 1:4)
hmc
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
