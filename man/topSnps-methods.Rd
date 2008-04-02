\name{topSnps-methods}
\docType{methods}
\alias{topSnps-methods}
\alias{topSnps}
\alias{topSnps,cwSnpScreenResult-method}
\title{ report on the smallest p-values for SNP-expression associations}
\description{ report on the smallest p-values for SNP-expression associations}
\section{Methods}{
\describe{

\item{x = "cwSnpScreenResult"}{ return a data frame with the
smallest p-values.  An optional argument \code{which = 'p.2df'} allows
selection of 2-d.f. test results.}
}}
\examples{
 data(hmceuB36.2021)
 g1 = gwSnpScreen(genesym("CPNE1"), hmceuB36.2021, chrnum(20))
 topSnps(g1)
}
\keyword{methods}
