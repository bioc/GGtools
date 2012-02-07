\name{transManager-class}
\Rdversion{1.1}
\docType{class}
\alias{transManager-class}
\alias{show,transManager-method}
\alias{transTab,transManager-method}
\alias{topGenes}
\alias{topScores}
\alias{nthScores}
\alias{geneIndcol}
\alias{geneNames}

\title{Class \code{"transManager"}}
\description{
simple container for manager of transScores output
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("transManager", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{base}:}{Object of class \code{"list"} includes
ff references for scores and indices of genes corresponding
to scores, and other metadata about the run }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "transManager")}: simple reporter }
	 }
}
%\references{
%%%  ~~put references to the literature/web site here~~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("transManager")
}
\keyword{classes}
