\name{eqBox}
\alias{eqBox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
descriptive plot of expression against genotype for cisAssoc
results
}
\description{
descriptive plot of expression against genotype for cisAssoc
results
}
\usage{
eqBox(gene, snp, se, tf, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{
identifier of gene in SummarizedExperiment se, must be present
in rownames(se)
}
  \item{snp}{
identifier of variant in VCF referenced by tf
}
  \item{se}{
SummarizedExperiment instance
}
  \item{tf}{
TabixFile reference for a VCF file that has been bgzipped
and tabix-indexed
}
  \item{\dots}{
will embellish plot; xlab and ylab already taken care of
}
}
\details{
The chromosome names in the VCF and the seqlevelsStyle of the Summarized
Experiment must match
}
%\value{
%%%  ~Describe the value returned
%%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  # forthcoming
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
