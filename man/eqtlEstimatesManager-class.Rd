\name{eqtlEstimatesManager-class}
\Rdversion{1.1}
\docType{class}
\alias{eqtlEstimatesManager-class}
\alias{[,eqtlEstimatesManager,ANY,ANY,ANY-method}

\title{Class \code{"eqtlEstimatesManager"}}
\description{
management of out-of-memory (ff-based) resources for
eQTL association estimates and their standard errors
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("eqtlEstimatesManager", ...)}.
Objects include metadata and ff-based reference elements.
}
\section{Slots}{
  \describe{
    \item{\code{fflist}:}{Object of class \code{"list"}; each element
is an ff object for accessing disk-resident matrices; typically
one element per chromosome analyzed  }
    \item{\code{call}:}{Object of class \code{"call"} for auditing
the call that generated the instance }
    \item{\code{sess}:}{Object of class \code{"ANY"} sessionInfo for
auditing}
    \item{\code{exdate}:}{Object of class \code{"ANY"} creation date }
    \item{\code{shortfac}:}{Object of class \code{"numeric"} multiplier
to allow representation of fractional components in short int }
    \item{\code{geneanno}:}{Object of class \code{"character"} name
of package that can be used to resolve probe identifiers }
    \item{\code{df}:}{Object of class \code{"numeric"} number of
degrees of freedom of each eQTL test }
    \item{\code{summaryList}:}{Object of class \code{"list"} \code{\link[snpStats]{col.summary}} is applied for each chromosome, so we can retrieve MAF
here }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eqtlTestsManager}"}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "eqtlEstimatesManager", i = "ANY", j = "ANY", drop = "ANY")}: \code{drop} is used to pick the 3rd dimension of the
array, is 1 for estimates, 2 for s.e. }
	 }
}
%\references{
%%%  ~~put references to the literature/web site here~~
%}
%%\author{
%%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("eqtlEstimatesManager")
}
\keyword{classes}
