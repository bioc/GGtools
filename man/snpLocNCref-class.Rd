\name{snpLocNCref-class}
\docType{class}
\alias{snpLocNCref-class}
\alias{getSnpData,snpLocNCref-method}
\alias{show,snpLocNCref-method}

\title{Class "snpLocNCref" container for netCDF reference
 to snp location data}
\description{Class "snpLocNCref" container for netCDF reference
 to snp location data}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("snpLocNCref", ...)}.
 They extend class character.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"character"} intended
to be the pathname of the netCDF store for snp location information.
This store will have two variables: chr (integer denoting
chromosome, in humans 23=X, 24=Y), and cumloc (numeric
cumulative base pair count to SNP location). }
  }
}
\section{Extends}{
Class \code{"\linkS4class{character}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "character", distance 2.
Class \code{"\link[Biobase:class.characterORmiame]{characterORMIAME}"}, by class "character", distance 2.
}
\section{Methods}{
  \describe{
    \item{getSnpData}{\code{signature(x = "snpLocNCref")}: returns
list with elements \code{chr} and \code{cumloc} }
    \item{show}{\code{signature(object = "snpLocNCref")}: returns
concise string as pathname may be long. }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("snpLocNCref")
}
\keyword{classes}
