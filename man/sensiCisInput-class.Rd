\name{sensiCisInput-class}
\Rdversion{1.1}
\docType{class}
\alias{sensiCisInput-class}
\alias{sensanal,sensiCisInput,numeric-method}
\alias{show,sensiCisInput-method}

\title{Class \code{"sensiCisInput"}}
\description{
Manage references to collections of cis-eQTL searches
for sensitivity analysis.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sensiCisInput", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{cisMgrFiles}:}{Object of class \code{"character"}:
a vector of filenames, each file is an instance of class
\code{\link{mcwBestCis-class}}  }
    \item{\code{cisMgrProperties}:}{Object of class \code{"list"} 
one vector with named elements per element of \code{cisMgrFiles}, with
components \code{rad}, \code{excl}, \code{maf}, \code{nperm},
\code{npc}; see details below.  }
    \item{\code{probeannopk}:}{Object of class \code{"character"},
identifying a bioconductor probe annotation package that can be used
to map probe identifiers to other vocabularies or
feature value sets   }
  }
}
\section{Methods}{
  \describe{
    \item{sensanal}{\code{signature(object = "sensiCisInput", fdrbound = "numeric")}: generates an instance of \code{\link{sensiCisOutput-class}} with
summarization of sensitivities }
    \item{show}{\code{signature(object = "sensiCisInput")}: concise
rendering }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%  ~~who you are~~
%}
\note{
This version of sensitivity analysis support is rudimentary
and involves manual construction of metadata that should
be extractable from analysis outputs.  The radius
of the cis search (and radius of excluded interior if
used) are identified as elements named \code{rad}
and \code{excl} in the \code{cisMgrProperties} vectors;
additional elements \code{maf}, \code{nperm}, and \code{npc}
define the lower bound for minor allele frequency, number of
permutations for plug-in FDR computation, and number of
principal components removed to adjust for expression heterogeneity
in the associated cis-eQTL search.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("sensiCisInput")
}
\keyword{classes}
