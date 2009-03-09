\name{aafSNP-class}
\Rdversion{1.1}
\docType{class}
\alias{aafSNP-class}
\alias{aafSNP}

\title{Class "aafSNP" -- container for HTML rendering of SNP metadata }
\description{Class "aafSNP" -- container for HTML rendering of SNP metadata }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("aafSNP", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"character"} will typically hold rs ids from dbSNP }
  }
}
\section{Extends}{
Class \code{"\linkS4class{character}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "character", distance 2.
Class \code{"\linkS4class{data.frameRowLabels}"}, by class "character", distance 2.
Class \code{"\linkS4class{characterORMIAME}"}, by class "character", distance 2.
Class \code{"\linkS4class{ANYTHING}"}, by class "character", distance 2.
Class \code{"\linkS4class{characterORNULL}"}, by class "character", distance 2.
Class \code{"\linkS4class{vectorORfactor}"}, by class "character", distance 3.
}
\section{Methods}{
The constructor has the same name, and operates on a list of character vectors.
It is expected that you would have a vector of rs numbers  for each gene,
thus a list of vectors with elements corresponding to genes.

getURL, getHTML are defined; see \code{\link[annaffy]{getURL}}, for example.
Apes the handling of UniGene links.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("aafSNP")
}
\keyword{classes}
