\name{maxchisq-class}
\Rdversion{1.1}
\docType{class}
\alias{maxchisq-class}
\alias{maxchisq}
\alias{min_p_vals}
\alias{min_p_vals,maxchisq,character,character,numeric-method}
\alias{show,maxchisq-method}

\title{Class "maxchisq" }
\description{container for results of cis-trans eQTL searches, and a p-value extractor}	 
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("maxchisq", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} currently representation
is simple -- a named list of named vectors of chisquared statistics corresponding to SNP,
a value for the d.f. of the chisq stats, the gene for which chisq was maximized for each
SNP, and some production metadata.  Note that a type parameter allows computation of 
max chisq stats per SNP (over genes) or per gene (over SNP)}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
Class \code{"\link[IRanges]{vectorORfactor-class}"}, by class "list", distance 3.
}
\section{Methods}{
  \describe{
    \item{min_p_vals}{\code{signature(mcs = "maxchisq", mtcorr = "character", type = "character", sidedness="numeric")}: 
\code{mtcorr} is the \code{proc} token for \code{\link[multtest]{mt.rawp2adjp}}.   Specifically,
if mtcorr is set to \code{"BH"}, the Benjamini-Hochberg FDR transformation is applied.
If mtcorr is set to \code{"none"}, nothing is done.

\code{type} determines the scope of the corrections.  Options are \code{""} which must be
used if \code{mtcorr} is \code{"none"}, \code{"chr_specific"}, with which the testing
corrections are made within chromosomes, or \code{"global"}, with which
the testing corrections are made over all tests over the whole genome.

sidedness determines whether a 2 sided (2*(1-pchisq)) or 1 sided
p-value is returned.  supply the factor 1 or 2 as desired.

 }
    \item{show}{\code{signature(object = "maxchisq")}: concise but informative report }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("maxchisq")
# also see example(diagffCC) for illustrations
}
\keyword{classes}
