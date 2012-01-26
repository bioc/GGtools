\name{geneIndcol}
\alias{geneIndcol}
\alias{geneNames}
\alias{locusNames}
\alias{nthScores}
\alias{topGenes}
\alias{topScores}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
tools for working with transManager instances
}
\description{
tools for working with transManager instances
}
\usage{
geneIndcol(tm, col)
geneNames(tm)
locusNames(tm)
nthScores(tm, n)
topGenes(tm)
topScores(tm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tm}{
instance of \code{\link{transManager-class}}
}
  \item{col}{
column selector
}
  \item{n}{
column selector
}
}
\details{
transManager instances have two ff matrices
of size L x K where L is the number of SNP
and K is the number of best feature scores
to be retained.  One matrix holds the scores,
the other holds the indices of the gene list
identifying the genes yielding the associated
scores.  Rows of the scores matrix
are sorted; the leftmost column of the scores
matrix is the maximum score.
}
\value{
each function returns a vector
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tm, col) 
tm@base$inds[, col]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
