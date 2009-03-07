\name{cisSnpTests}
\Rdversion{1.1}
\alias{cisSnpTests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
perform tests for eQTL cis to specified genes
}
\description{
perform tests for eQTL cis to specified genes
}
\usage{
cisSnpTests(fmla, smls, radius, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fmla}{
standard formula. LHS can be a GeneSet with AnnotationIdentifier
geneIdType.  RHS can be predictor formula component using variables
in pData of \code{smls}
}
  \item{smls}{
instance of smlSet
}
  \item{radius}{
numeric value: number of bases up and downstream from probe
CHRLOC to be examined for SNP
}
  \item{\dots}{ not in use
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a list of cwSnpScreen instances
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
Getting SNP locations is slow for the first event while metadata
are brought into scope.  Subsequent calls are faster.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(GSEABase)
# two genes on chr 20
gs1 = GeneSet(c("CPNE1", "ADA"), geneIdType = SymbolIdentifier())
gs2 = gs1
organism(gs2) = "Homo sapiens"
geneIdType(gs2) = AnnotationIdentifier("illuminaHumanv1.db")
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
cc = cisSnpTests(gs2~male, hmceuB36.2021, radius=1e5)
lapply(cc, function(x) length(p.value(x@.Data[[1]])))
cc = cisSnpTests(gs2~male, hmceuB36.2021, radius=1e6)
lapply(cc, function(x) length(p.value(x@.Data[[1]])))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
