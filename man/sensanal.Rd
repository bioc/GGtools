\name{sensanal}
\alias{sensanal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize information from a collection of eQTL searches
for sensitivity assessment
}
\description{
Summarize information from a collection of eQTL searches
for sensitivity assessment
}
\usage{
sensanal(object, fdrbound)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
instance of \code{\link{sensiCisInput-class}}
}
  \item{fdrbound}{
numeric upper bound on FDR for declarations of eQTL yield
%%     ~~Describe \code{fdrbound} here~~
}
}
\details{
Sensitivity analysis for cis-eQTL search involves checking
effects of scope of search, allele frequency filtering, and
adjustment for expression heterogeneity on eQTL declarations.
In this version, we focus on collections of outputs of
\code{\link{best.cis.eQTLs}}, to which the values
of tuning parameters are bound.  These collections
are identified in a \code{\link{sensiCisInput-class}}
instance, and the \code{sensanal} function processes
these outputs into a \code{\link{sensiCisOutput-class}} instance
for tabulation and visualization.
}

\value{
a \code{\link{sensiCisOutput-class}} instance
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
