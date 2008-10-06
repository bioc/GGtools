\name{hbTestResults-class}
\docType{class}
\alias{hbTestResults-class}
\alias{pvals,hbTestResults-method}

\title{Class "hbTestResults" holds results of tests
  of association of expression levels with
  haplotype within haplotype block}
\description{Class "hbTestResults" holds results of tests
  of association of expression levels with
  haplotype within haplotype block}
   
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("hbTestResults", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{hscores}:}{Object of class \code{"list"} series of
haplo.stats:::haplo.score results for blocks  }
    \item{\code{locs}:}{Object of class \code{"numeric"} locations
at which blocks were found (mean location within each block) }
    \item{\code{chrnum}:}{Object of class \code{"chrnum"} chromosome
being analyzed }
    \item{\code{smlSetName}:}{Object of class \code{"character"} name
of the smlSet harboring data in use }
    \item{\code{rsid}:}{Object of class \code{"ANY"} can be
a dbSNP id to use as an anchor, or a number constituting
absolute chromosomal location at which blocks will be
sought }
    \item{\code{rad}:}{Object of class \code{"numeric"} radius
in base pairs around the \code{rsid} to be searched for blocks }
    \item{\code{ldStruc}:}{Object of class \code{"ANY"} the result
of the mapLD:::mapLD function }
  }
}
\section{Methods}{
  \describe{
    \item{pvals}{\code{signature(x = "hbTestResults")}: extracts
p-values for global score tests, one per block }
    \item{locs}{\code{signature(x = "hbTestResults")}: extracts
locations of haplotype blocks found  (average SNP location within block) }
    \item{hscores}{\code{signature(x = "hbTestResults")}: extracts
\code{\link{haplo.score}} results as a list, for all blocks}
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu>}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("hbTestResults")
}
\keyword{classes}
