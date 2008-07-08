\name{gwSnpScreenResult-class}
\docType{class}
\alias{gwSnpScreenResult-class}
\alias{multiGwSnpScreenResult-class}
\alias{cwSnpScreenResult-class}
\alias{plot,gwSnpScreenResult,ANY-method}
\alias{plot,cwSnpScreenResult,ANY-method}
\alias{plot,gwSnpScreenResult-method}
\alias{plot,cwSnpScreenResult-method}
\alias{show,gwSnpScreenResult-method}
\alias{getAbsSnpLocs,cwSnpScreenResult-method}
\alias{getAbsSnpLocs}
\alias{show,multiGwSnpScreenResult-method}
\alias{[,gwSnpScreenResult-method}
\alias{[,gwSnpScreenResult,ANY,ANY-method}
\alias{[,gwSnpScreenResult,ANY,ANY,ANY-method}
\alias{getSnpOffsets,cwSnpScreenResult-method}
\alias{getSnpOffsets}


\title{Class "gwSnpScreenResult" -- containers for
genome-wide and chromosome-wide snp screen statistics }
\description{Class "gwSnpScreenResult" -- containers for
genome-wide and chromosome-wide snp screen statistics }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gwSnpScreenResult", ...)}.
Mild extensions to 'list'.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} list of
data.frames or matrices of test results }
    \item{\code{gene}:}{Object of class \code{"character"} gene name }
    \item{\code{psid}:}{Object of class \code{"character"} probe identifier }
    \item{\code{annotation}:}{Object of class \code{"character"} -- propagated
from smlSet }
    \item{\code{snpLocExtRef}:}{Object of class \code{"character"}, full pathname
of external store for snp locations relevant to the smlSet instance from
which gwSnpScreenResult was computed.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "gwSnpScreenResult", y = "ANY")}: 
plots minus log p values of association tests along chromosomal
coordinates}
    \item{show}{\code{signature(object = "gwSnpScreenResult")}: brief report }
    \item{[}{\code{signature(object = "gwSnpScreenResult")}: restricts
chromosomes represented }
    \item{getSnpOffsets}{\code{signature(object = "cwSnpScreenResult")}: 
     return vector of offsets from telomere end to first SNP per
chromosome }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("gwSnpScreenResult")
data(hmceuB36.2021)
objects("package:GGBase")
x = gwSnpScreen(genesym("CPNE1"), hmceuB36.2021, chrnum(20))
x
}
\keyword{classes}
