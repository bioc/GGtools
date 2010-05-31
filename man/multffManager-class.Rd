\name{multffManager-class}
\Rdversion{1.1}
\docType{class}
\alias{multffManager-class}
\alias{show,multffManager-method}
\alias{coerce,multffManager,eqtlTestsManager-method}
\alias{[,multffManager,missing,probeId,ANY-method}
\alias{[,multffManager,rsid,probeId,ANY-method}
\alias{[,multffManager,rsid,missing,ANY-method}

\title{Class "multffManager"  }
\description{	 coordinates access to and interrogation of multipopulation eQTL searches }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("multffManager", ...)}.
These extend list during the experimental development phase.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
Class \code{"\linkS4class{vectorORfactor}"}, by class "list", distance 3.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "multffManager")}: concise report that
provides an excerpt from the ff image }
    \item{[}{\code{signature(x = "multffManager"), i, j, ...}: you can extract results by rsid or probeId with customary
bracket semantics, with the exception that if the SNP request spans multiple chromosomes, you will get a list of results }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu>}
\note{ 
> names(dd)
%[1] "fflist"       "call"         "runname"      "targdir"      "generangetag"
%[6] "filenames"    "df"           "vmode"        "shortfac"    
Currently components of .Data are 
\describe{
 \item{\code{fflist}}{a list of ff references, to tables holding sums of chi-squared statistics
accumulated across populations}
 \item{\code{call}}{for auditing, the initial call}
 \item{\code{runname}}{an arbitrary user-supplied tag}
 \item{\code{targdir}}{the folder used to write the ff files}
 \item{\code{generangetag}}{a generated tag giving the scope of the gene set used for searches}
 \item{\code{filenames}}{a character vector of the ff file paths}
 \item{\code{df}}{numeric value of the number of populations summed}
 \item{\code{vmode}}{ff specification of virtual mode of data values; if 'short', rescale using
shortfac}
 \item{\code{shortfac}}{factor by which chisquared deviates were multiplied so that a short
int can represent without too much coarsening}
}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
#
# seems to throw file error in CMD check on windows
#
if (.Platform$OS.type != "windows") {
  example("multffCT")
  dem
  getClass(class(dem))
  dem$fflist[[1]]
  dem$df
  dem$filenames
  dem$vmode
  dem$call
  }
}
\keyword{classes}
