\name{eqtlTestsManager-class}
\Rdversion{1.1}
\docType{class}
\alias{eqtlTestsManager-class}
\alias{[,eqtlTestsManager,ANY,ANY,ANY-method}
\alias{show,eqtlTestsManager-method}
\alias{topFeats,probeId,eqtlTestsManager-method}
\alias{topFeats,rsid,eqtlTestsManager-method}

\title{Class \code{"eqtlTestsManager"}}
\description{
manage out-of-memory elements of an eQTL search
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("eqtlTestsManager", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fffile}:}{Object of class \code{"ff_matrix"} chisquared
statistics stored as short ints in ff out of memory file }
    \item{\code{call}:}{Object of class \code{"call"} audit
of creation call }
    \item{\code{sess}:}{Object of class \code{"ANY"} session info
structure at time of creation}
    \item{\code{exdate}:}{Object of class \code{"ANY"} date
at time of creation }
    \item{\code{shortfac}:}{Object of class \code{"numeric"} number
by which chisq stats are multiplied to allow recovery of precision }
    \item{\code{geneanno}:}{Object of class \code{"character"} string
naming annotation package relevant for probe identifier translation }
    \item{\code{df}:}{Object of class \code{"numeric"} degrees
of freedom of chisq stats }
    \item{\code{summaryList}:}{Object of class \code{"list"} 
list of genotype statistical summaries }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "eqtlTestsManager", i = "ANY", j = "ANY", drop = "ANY")}: extract chisq statistics properly rescaled from short int
to double }
    \item{show}{\code{signature(object = "eqtlTestsManager")}: concise
report }
    \item{topFeats}{\code{signature(feat = "probeId", mgr = "eqtlTestsManager")}: extract highest scores for SNP associated with given probeId }
    \item{topFeats}{\code{signature(feat = "rsid", mgr = "eqtlTestsManager")}: extract highest scores for probes associated with given SNP }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
instances are created by \code{\link{eqtlTests}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("eqtlTestsManager")
}
\keyword{classes}
