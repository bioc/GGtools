\name{transTab}
\alias{transTab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
tabulate results of transScores run
}
\description{
tabulate results of transScores run
}
\usage{
transTab(x, snps2keep=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a transManager instance.
}
  \item{snps2keep}{
character vector used for filtering snps whose scores will be retained;
intersection with snps named in x will be used.
}
%  \item{snpchr}{
%%%     ~~Describe \code{snpchr} here~~
%string denoting the chromosome from which SNP genotypes were drawn for trans-association with gene expression
%}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
data.frame instance
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
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
