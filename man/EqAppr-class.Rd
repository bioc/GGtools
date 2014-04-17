\name{EqAppr-class}
\Rdversion{1.1}
\docType{class}
\alias{EqAppr-class}
\alias{calfig,EqAppr,character-method}
\alias{getModnames,EqAppr-method}
\alias{getPruned,EqAppr-method}
\alias{getSens,EqAppr-method}
\alias{getUnpruned,EqAppr-method}
\alias{show,EqAppr-method}

\title{Class \code{"EqAppr"}}
\description{
Manage the appraisal of an eQTL search
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EqAppr", ...)}, or via \code{buildEqAppr()}
}
\section{Slots}{
  \describe{
    \item{\code{meta}:}{Object of class \code{"ApprMeta"} basic descriptive information
about source analysis }
    \item{\code{sens}:}{Object of class \code{"ApprSens"} outputs of sensitivity analysis }
    \item{\code{pruned}:}{Object of class \code{"ApprRes"} outputs of appraisal for LD-pruned predictions }
    \item{\code{unpruned}:}{Object of class \code{"ApprRes"} outputs of general appraisal without LD pruning }
  }
}
\section{Methods}{
  \describe{
    \item{calfig}{\code{signature(x = "EqAppr", ind = "character")}: ... }
    \item{getModnames}{\code{signature(x = "EqAppr")}: ... }
    \item{getPruned}{\code{signature(x = "EqAppr")}: ... }
    \item{getSens}{\code{signature(x = "EqAppr")}: ... }
    \item{getUnpruned}{\code{signature(x = "EqAppr")}: ... }
    \item{show}{\code{signature(object = "EqAppr")}: ... }
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
showClass("EqAppr")
}
\keyword{classes}
