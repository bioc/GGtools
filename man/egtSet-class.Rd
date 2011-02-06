\name{egtSet-class}
\Rdversion{1.1}
\docType{class}
\alias{egtSet-class}
\alias{getSNP,egtSet,ANY,ScalarCharacter-method}
\alias{show,egtSet-method}
\alias{convertSS}

\title{Class "egtSet"}
\description{
expression plus genotype set, with external representation of
genotypes
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("egtSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.  A convenience function
convertSS is provided and described below.
}
\section{Slots}{
  \describe{
    \item{\code{gtlist}:}{Object of class \code{"list"} a list of
ff instances, intended to have one element per chromosome }
    \item{\code{snpNames}:}{Object of class \code{"list"} list
of strings, one element per chromosome; each element has length
the number of SNP on the associated chromosome }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} eSet
container for expression values }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} 
as in eSet }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} 
as in eSet }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} as in eSet }
    \item{\code{annotation}:}{Object of class \code{"character"} as in eSet }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} as in eSet }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} as in eSet }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{getSNP}{\code{signature(egts = "egtSet", chrind = "ANY", tok = "ScalarCharacter")}: retrieves matrix of raw corresponding to genotypes for
one SNP over all samples }
    \item{show}{\code{signature(object = "egtSet")}: print method }
	 }
}
%\references{
%%%  ~~put references to the literature/web site here~~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
the \code{convertSS} function will create an egtSet instance from an
smlSet instance.  Arguments are smlSet, ffstub, the latter a string
that will be used to prefix names of per-chromosome (really per element
of  smList(smlSet)) ff files corresponding to SnpMatrix instances
}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("egtSet")
curd = getwd()
td = tempdir()
setwd(td)
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
fceu = convertSS(hmceuB36.2021)
fceu
dir()
fceu@gtlist
setwd(curd)
}
\keyword{classes}
