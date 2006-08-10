\name{ggExprSet-class}
\docType{class}
\alias{ggExprSet-class}
\alias{ggExprSet}
\alias{snps,ggExprSet-method}

\title{Class "ggExprSet", representation of genetical
genomics (expression + SNP) data }
\description{	 
Class "ggExprSet", representation of genetical
genomics (expression + SNP) data 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ggExprSet", ..., versions)}.
These are identical to \code{\link[Biobase]{exprSet-class}} instances
from Biobase, but additional methods are available.
}
\section{Slots}{
	 \describe{
    \item{\code{exprs}:}{Object of class \code{"exprMatrix"} as in \code{\link[Biobase]{exprSet-class}} }
    \item{\code{se.exprs}:}{Object of class \code{"exprMatrix"} as in \code{\link[Biobase]{exprSet-class}} }
    \item{\code{description}:}{Object of class \code{"characterORMIAME"} as in \code{\link[Biobase]{exprSet-class}} }
    \item{\code{annotation}:}{Object of class \code{"character"} as in \code{\link[Biobase]{exprSet-class}} }
    \item{\code{notes}:}{Object of class \code{"character"} as in \code{\link[Biobase]{exprSet-class}} }
    \item{\code{reporterInfo}:}{Object of class \code{"data.frameOrNULL"} as in \code{\link[Biobase]{exprSet-class}} }
    \item{\code{phenoData}:}{Object of class \code{"phenoData"} as in \code{\link[Biobase]{exprSet-class}}, but note that in this class, we use phenoData to encode genotype! }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} as in \code{\link[Biobase]{exprSet-class}} }
  }
}
\section{Extends}{
Class \code{\link[Biobase]{exprSet-class}}, directly.
Class \code{"annotatedDataset"}, by class \code{"exprSet"}.
Class \code{"VersionedBiobase"}, by class \code{"exprSet"}.
Class \code{"Versioned"}, by class \code{"exprSet"}.
}
\section{Methods}{
  \describe{
    \item{snps}{\code{signature(x = "ggExprSet")}: accessor to
the data.frame holding the genotypes }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }

\seealso{
	\code{\link[Biobase]{exprSet-class}} 
}
\examples{
getClass("ggExprSet")
# be patient...
data(c20GGceu)
snps(c20GGceu)[1:4,1:4]
}
\keyword{classes}
