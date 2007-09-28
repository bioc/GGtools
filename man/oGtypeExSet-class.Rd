\name{oGtypeExSet-class}
\docType{class}
\alias{oGtypeExSet-class}
\alias{coerce,oGtypeExSet,ExpressionSet-method}
\alias{initialize,oGtypeExSet-method}
\alias{show,oGtypeExSet-method}
\alias{addAnno}
\alias{getSNPindices}
\alias{tx_rsnum}
\title{Class "oGtypeExSet" ~~~ }
\description{	 Experimental class for combining oligo
SnpCallSet with ExpressionSet }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("oGtypeExSet", 
snpCalls=list(), dbConns=list(), x)}
where \code{x} is an \code{\link[Biobase]{ExpressionSet}} instance.
}
\section{Slots}{
	 \describe{
    \item{\code{snpCalls}:}{Object of class \code{"list"} named
list of \code{\link[oligo]{SnpCallSet}} instances  }
    \item{\code{dbConns}:}{Object of class \code{"list"} list
of SQLite connections for SNP annotation structures }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} from expression data }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{\link[Biobase]{ExpressionSet-class}}, directly.
Class \code{\link[Biobase]{eSet-class}}, by class "ExpressionSet", distance 2.
Class \code{\link[Biobase]{VersionedBiobase-class}}, by class "ExpressionSet", distance 3.
Class \code{\link[Biobase]{Versioned-class}}, by class "ExpressionSet", distance 4.
}
\section{Methods}{
No methods defined with class "oGtypeExSet" in the signature.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\note{ addAnno, getSNPindices, and tx.rsnum need further
development.}

% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%}
\examples{
# use new("oGtypeExSet", snpCalls=list(Nsp=nspCalls, Sty=styCalls),
#    dbConns=list(), x3pSet)
showClass("oGtypeExSet")
}
\keyword{classes}
