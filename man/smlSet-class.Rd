\name{smlSet-class}
\docType{class}
\alias{hmceuB36.2021}
\alias{smlSet-class}
\alias{chrnum-class}
\alias{rsNum-class}
\alias{[,smlSet-method}
\alias{exprs,smlSet-method}
%\alias{gwSnpScreen,genesym,smlSet-method}
\alias{nsamp,smlSet-method}
\alias{rawSNP,smlSet,numeric-method}
\alias{sampleNames,smlSet-method}
\alias{show,smlSet-method}
\alias{smList,smlSet-method}
\alias{smlEnv,smlSet-method}
\alias{getSnpLocs,smlSet-method}
\alias{getSnpLocs}
\alias{getSnpChroms,smlSet-method}
\alias{getSnpChroms}
\alias{getSnpData}
\alias{getSnpData,character-method}
\alias{rawSNP}
\alias{smList}
\alias{chrnum}
\alias{rsNum}
\alias{chrnum,numeric-method}
\alias{rsNum,character-method}
\alias{getSnpChroms,smlSet,logical-method}
\alias{getSnpChroms,smlSet,missing-method}
\alias{getSnpData,character,character-method}
\alias{getSnpLocs,smlSet,logical-method}
\alias{getSnpLocs,smlSet,missing-method}
%\alias{snps,racExSet,missing-method}
%\alias{snps,smlSet,chrnum-method}

\title{Class "smlSet" -- efficient representation of SNP+expression data,
based on Clayton snpMatrix package}
\description{Efficient representation of SNP+expression data,
based on Clayton snpMatrix package.  Supports efficient subsetting
of whole-genome snp genotyping by chromosome.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("smlSet", smlEnv, assayData, phenoData, featureData, experimentData, annotation, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{smlEnv}:}{Object of class \code{"environment"} -- there is
one key 'smList' with value a list of \code{snp.matrix} from snpMatrix package }
    \item{\code{snpLocPathMaker}:}{Object of class \code{"function"} --
a function of one argument that returns a full pathname to the netCDF
store for SNP location information }
    \item{\code{chromInds}:}{Object of class \code{"numeric"}, a numeric
vector telling what chromosomes are represented. If chromosomes 20 and 21
are present, the values is c(20,21).  X and Y are 23, 24.}
    \item{\code{organism}:}{Object of class \code{"character"}; string such as Hs for
modern humans. }
    \item{\code{assayData}:}{Object of class \code{"AssayData"}; container for
expression data as in ExpressionSet }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}; container
for sample data as in ExpressionSet  }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}; container
for expression feature data as in ExpressionSet}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}; container for MIAME
schema  }
    \item{\code{annotation}:}{Object of class \code{"character"}; vector with
names \code{exprs} and \code{snps}; \code{exprs} should be the name of a
SQLite-based Bioconductor annotation package  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "smlSet")}: subsetting must be done with care.
The usual idiom is \code{X[G, S]} where \code{G} restricts reporters and
\code{S} restricts samples.  Reporter restriction:
At present one can use an instance of \code{chrnum} class to obtain
a subsetted object with SNP data restricted to specified chromosomes.
Sample restriction: works in the customary way for \code{eSet} instances. }
    \item{exprs}{\code{signature(object = "smlSet")}: return the
expression data matrix. }
    \item{gwSnpScreen}{\code{signature(sym = "genesym", sms = "smlSet")}: Execute
snpMatrix snp.rhs.tests for all SNP in the genome.  Returns instance
of gwSnpScreenResult. }
    \item{cwSnpScreen}{\code{signature(sym = "genesym", sms = "smlSet",
      cnum="chrnum")}: Execute
snpMatrix snp.rhs.tests for all SNP on the specific
 chromosome.  Returns instance
of cwSnpScreenResult. }
    \item{nsamp}{\code{signature(x = "smlSet")}: returns number of samples. }
    \item{rawSNP}{\code{signature(x = "smlSet", chrind = "numeric")}: returns
byte (raw) representation of SNP data in snp.matrix format for specified
chromosome.  }
    \item{getSnpLocs}{\code{signature(x = "smlSet")}: returns
all relevant locations in order of rs numbers in the
snp.matrices as a single vector}
    \item{getSnpChroms}{\code{signature(x = "smlSet")}: returns
all relevant chromosome assignments in order of rs numbers in the
snp.matrices as a single vector}
    \item{sampleNames}{\code{signature(object = "smlSet")}: returns
vector of sample names. }
    \item{show}{\code{signature(object = "smlSet")}: prints short report. }
    \item{smList}{\code{signature(x = "smlSet")}: returns the list of
snp.matrix encoding genotypes. }
    \item{smlEnv}{\code{signature(x = "smlSet")}: returns the environment
in which the smList list is contained. }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu>}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
#showClass("smlSet")
data(hmceuB36.2021)
rawSNP(hmceuB36.2021, 20)[1:5,1:5]
exprs(hmceuB36.2021)[1:5,1:5]
}
\keyword{classes}
