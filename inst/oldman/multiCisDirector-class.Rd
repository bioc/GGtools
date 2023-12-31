\name{multiCisDirector-class}
\Rdversion{1.1}
\docType{class}
\alias{multiCisDirector-class}
\alias{makeDiagDirector}
\alias{show,multiCisDirector-method}

\title{Class "multiCisDirector"}
\description{
manage multiple eqtlTestsManager instances, typically as interim results from a 
run of cisProxScores
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("multiCisDirector", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{mgrs}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "multiCisDirector")}: ... }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
\code{makeDiagDirector} is a tool that will generate all same-chromosome
eqtlTests from an smlSet instance or package and will create a director of this type.

makeDiagDirector and this class will be removed in GGtools 4.x;
use eqtlTests and best.cis.eQTLs for comparable functionality.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cisProxScores}}
}
\examples{
showClass("multiCisDirector")
}
\keyword{classes}
