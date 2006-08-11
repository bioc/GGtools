
\name{racExSet-class}
\docType{class}
\alias{racExSet-class}
\alias{snps,racExSet-method}
\alias{snpNames,racExSet-method}
\alias{racAssays,racExSet-method}
\alias{litsn}
\alias{litex}
\alias{pdd}
\alias{racAssays}
\alias{snpNames}


\title{Class "racExSet" for combining RareAlleleCount representations
of SNPs, gene expression data, and other phenotype data }
\description{	 A container for combining RareAlleleCount representations 
of SNPs, gene expression data, and other phenotype data }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("racExSet", phenoData, experimentData, annotation, exprs, racs, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{racAssays}:}{Object of class \code{"AssayData"} SxN matrix
of 0, 1, 2, or NA, rownames are RS numbers from dbSNP, column names
are sampleNames }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} GxN matrix
of expression values }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} standard
phenoData instance }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} experiment
level metadata }
    \item{\code{annotation}:}{Object of class \code{"character"} tag }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} R class version
control metadata}
  }
}
\section{Extends}{
Class \code{"eSet"}, directly.
Class \code{"VersionedBiobase"}, by class \code{"eSet"}.
Class \code{"Versioned"}, by class \code{"eSet"}.
}
\section{Methods}{
  \describe{
    \item{snps}{\code{signature(x = "racExSet")}: accessor for the
racAssay matrix }
    \item{snpNames}{\code{signature(x = "racExSet")}: give RS numbers
for rows of racAssay matrix }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~
%
\seealso{
\code{\link{make_racExSet}} for construction; \code{\link[Biobase]{eSet}} ancestor
class
}
\examples{
# as in make_racExSet
data(litsn)
data(litex)
data(pdd)
dem = make_racExSet( exprs=litex, racs=t(litsn), pd=pdd, mi=new("MIAME"),
  anno="demo")
dem
snpNames(dem)[1:4]
featureNames(dem)[1:4]
}
\keyword{classes}
