\name{gwSnpTests}
\alias{gwSnpTests}
\alias{gwSnpTests,formula,smlSet,missing-method}
\alias{gwSnpTests,formula,smlSet-method}
\alias{topSnps}
\alias{topSnps,gwSnpScreenResult-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
execute a series of
tests for association between genotype and expression
}
\description{
execute a series of tests for association between genotype and expression
}
\usage{
gwSnpTests(sym, sms, ...)
topSnps(x, n=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sym}{
instance of \code{\linkS4class{probeId}} or
\code{\linkS4class{genesym}}
}
  \item{sms}{
instance of \code{\link[ggbase2]{smlSet-class}}
}
  \item{x}{instance of \code{gwSnpScreenResult}}
  \item{n}{integer, number of test results to be reported,
sorted decreasing from the most significant}
  \item{\dots}{
%arguments passed to \code{\link[snpStats]{snp.rhs.tests}}
not used
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
an instance of the gwSnpScreenResult class,
to be examined by \code{topSnps}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
The most basic application yields one d.f. chi-squared
statistics based on score tests.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
s20 = getSS("GGtools", "20")
t1 = gwSnpTests(genesym("CPNE1")~male, s20)
topSnps(t1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
