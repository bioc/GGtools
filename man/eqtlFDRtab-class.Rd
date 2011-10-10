\name{eqtlFDRtab-class}
\Rdversion{1.1}
\docType{class}
\alias{eqtlFDRtab-class}
\alias{show,eqtlFDRtab-method}
\alias{show,gwScores-method}
\alias{show,pwScores-method}

\title{Class \code{"eqtlFDRtab"}}
\description{
manage results of genewiseFDRtab and allied methods
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("eqtlFDRtab", ...)}.
%%  ~~ describe objects here ~~ 
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
%Class \code{"\linkS4class{vectorORfactor}"}, by class "list", distance 3.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "eqtlFDRtab")}: ... }
    \item{show}{\code{signature(object = "gwScores")}: utility container en route to eqtlFDRtab }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("eqtlFDRtab")
}
\keyword{classes}
