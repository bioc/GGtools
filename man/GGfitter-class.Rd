\name{GGfitter-class}
\docType{class}
\alias{GGfitter-class}
\alias{fastAGM}
\alias{fastHET}
\alias{fastAGMfitter}
\alias{fastHETfitter}
\alias{snpScreen,racExSet,snpMeta,genesym,formula,GGfitter,ANY-method}

\title{Class "GGfitter" wrapper for special functions to do
 whole genome screens }
\description{	 wraps functions that call C code }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GGfitter", ...)}.
They can be supplied as arguments to \code{\link{snpScreen}}.
}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{Object of class \code{"character"} giving a name }
    \item{\code{func}:}{Object of class \code{"function"} the actual
function that calls C to perform the screen }
  }
}
\section{Methods}{
  \describe{
    \item{snpScreen}{\code{signature(racExSet = "racExSet", snpMeta = "snpMeta", gene = "genesym", formTemplate = "formula", fitter = "GGfitter", gran = "ANY")}: primary usage of
this class }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu>}
\note{ Previously, deparse(substitute(fn)) was used to get information
on the fitter.  This is not reliable, and all fitting functions used
in this context will have to be reflective instances of this class.}



%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("GGfitter")
}
\keyword{classes}
