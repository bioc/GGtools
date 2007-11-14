\name{snpMeta-class}
\docType{class}
\alias{snpMeta-class}
\alias{snpMetaWhole-class}
\alias{snpMetaWhole}
\alias{pos}
\alias{pos,snpMeta-method}
\alias{wrapSNPmetaWh}
\alias{chr20meta}
\alias{[,snpMeta-method}
\alias{chromosome,snpMeta-method}
\alias{coerce,snpMeta,data.frame-method}
\alias{dim,snpMeta-method}
\alias{nrow,snpMeta-method}
\alias{show,snpMeta-method}

\title{Class "snpMeta" -- HapMap (or Wellcome INBREDS) -based metadata structures for SNPs }
\description{	 HapMap (or Wellcome INBREDS) -based metadata structures for SNPs }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("snpMeta", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{meta}:}{Object of class \code{"environment"} holds the table in
the \code{meta} element,
in the form of a \code{data.frame}, with rownames giving the SNP ids,
and columns \code{pos}, \code{strand}}
    \item{\code{chromosome}:}{Object of class \code{"character"} string
identifying chromosome on which SNPs are documented }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "snpMeta")}: select subset using rs number}
    \item{chromosome}{\code{signature(x = "snpMeta")}: obtain chromosome name }
    \item{pos}{\code{signature(x = "snpMeta")}: obtain named vector of positions }
    \item{coerce}{\code{signature(from = "snpMeta", to = "data.frame")}: return
        data.frame with all metadata }
    \item{dim}{\code{signature(x = "snpMeta")}: give dimensions }
    \item{nrow}{\code{signature(x = "snpMeta")}: tell number of entries }
    \item{show}{\code{signature(object = "snpMeta")}: give abbreviated report of contents
with exemplar records}
	 }
}
\references{ \url{www.hapmap.org} , \url{www.well.ox.ac.uk/mouse/INBREDS} }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\link{CLASSNAME-class}} for links to other classes

\examples{
data(chr20meta)
chr20meta["rs6060535",]
}
\keyword{classes}
