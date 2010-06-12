\name{eqtlTestsManager-class}
\Rdversion{1.1}
\docType{class}
\alias{eqtlTestsManager-class}
\alias{cisTransDirector-class}
\alias{[,eqtlTestsManager,rsid,probeId,ANY-method}
\alias{[,cisTransDirector,character,character,ANY-method}
\alias{[,cisTransDirector,character,missing,ANY-method}
\alias{[,cisTransDirector,missing,character,ANY-method}
\alias{show,eqtlTestsManager-method}
\alias{show,cisTransDirector-method}
\alias{probeNames,eqtlTestsManager-method}
\alias{probeNames,cisTransDirector-method}
\alias{probeNames}

\title{Class "eqtlTestsManager"}
\description{
interface to ff files that store results for
large numbers of eQTL tests
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("eqtlTestsManager", ...)},
or \code{new("cisTransDirector", ...)}.  The \code{mkCisTransDirector} function
should be used for the latter task.

A manager object collects metadata and reference information regarding
tests relating 
a single set of expression measures (gene-oriented) and a collection
of structural variants (snp-oriented).

A director object collects metadata and reference information for a
specified set of managers.
}
\section{Slots}{
  \describe{
    \item{\code{fflist}:}{Object of class \code{"list"} ~~ }
    \item{\code{call}:}{Object of class \code{"call"} ~~ }
    \item{\code{sess}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{exdate}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{shortfac}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{geneanno}:}{Object of class \code{"character"} ~~ }
    \item{\code{df}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "eqtlTestsManager", i = "rsid", j = "probeId", drop = "ANY")}: This
gives matrix-like extraction idiom to retrieve chisquared statistics from the ff
archives for eQTL searches }
    \item{[}{\code{signature(x = "cisTransDirector", i = "character", j = "character", drop = "ANY")}: ... }
    \item{show}{\code{signature(object = "eqtlTestsManager")}: ... }
    \item{show}{\code{signature(object = "cisTransDirector")}: ... }
    \item{probeNames}{\code{signature(object = "eqtlTestsManager")}: extract the probe names as
a vector }
    \item{probeNames}{\code{signature(object = "cisTransDirector")}: extract the probe names as a list
with one element per manager }
	 }
}
%\references{
%%%  ~~put references to the literature/web site here~~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
Instances of this class can be coerced to instances of eqtlTestsManager to facilitate management
by a cisTransDirector
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
# look at example(eqtlTests) for workout
showClass("eqtlTestsManager")
showClass("cisTransDirector")
}
\keyword{classes}
