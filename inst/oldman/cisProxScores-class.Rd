\name{cisProxScores-class}
\Rdversion{1.1}
\docType{class}
\alias{cisProxScores-class}
\alias{show,cisProxScores-method}

\title{Class "cisProxScores"}
\description{
extends list to manage collections of eQTL test scores
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cisProxScores", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
    \item{\code{call}:}{Object of class \code{"call"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
Class \code{vectorORfactor}, by class "list", distance 3.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "cisProxScores")}: concise report }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("cisProxScores")
}
\keyword{classes}
