
\name{racExSet-class}
\docType{class}
\alias{racExSet-class}
\alias{genesym-class}
\alias{genesym}
\alias{getpsid}
\alias{snps,racExSet-method}
\alias{snpNames,racExSet-method}
\alias{racAssays,racExSet-method}
\alias{oneFit,racExSet-method}
\alias{litsn}
\alias{oneFit}
\alias{litex}
\alias{pdd}
\alias{racAssays}
\alias{snpNames}
\alias{SNPalleles}
\alias{rarebase}
\alias{SNPalleles,racExSet-method}
\alias{rarebase,racExSet-method}
\alias{oneFit,racExSet,genesym,character,function-method}
\alias{oneFit,racExSet,character,character,function-method}
\alias{oneFit,racExSet,list,character,function-method}
\alias{oneFit,racExSet,list,formula,function-method}
\alias{show,racExSet-method}

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
    \item{\code{rarebase}:}{character vector of length G stating the
rare base of the SNP}
    \item{\code{SNPalleles}:}{character vector of length G stating the
alleles in the form "A/B"}
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
    \item{racAssays}{\code{signature(x = "racExSet")}: accessor for the
racAssays slot}
    \item{oneFit}{\code{signature(racExSet = "racExSet", geneid = "genesym", snpid = "character", fitfun = "function")}: use a genesym instance (e.g.,
\code{genesym("DDR1")} to identify the expression values,
an rs number to identify the snp genotype information, and
fitfun to fit a model with expression as dependent variable and
genotype (rare allele counts) as predictor.  Any fitting
function accepting a formula and data frame should succeed (e.g., \code{lm}). }
    \item{oneFit}{\code{signature(racExSet = "racExSet", geneid = "character", snpid = "character", fitfun = "function")}: as above, but can avoid the lookup of the
gene symbol by supplying the probe set identifier directly.  This is
problematic at present as most probe set identifier strings are not
valid R variable names, and so cannot serve directly in the formula.
Some backtick coding is needed to get this to work generally. }
    \item{oneFit}{\code{signature(racExSet = "racExSet", geneid = "list", snpid = "character", fitfun = "function")}: the \code{geneid} parameter must
be a named list of length one with content a numeric vector of
expression values conforming to the snp values.  Avoids lookup
of gene name and of expression values.}
   \item{oneFit}{\code{signature(racExSet = "racExSet", geneid = "list" , snpid="formula", fitfun="function")}: the \code{geneid} parameter must be bound to a named list of length 1
with content a numeric vector
with one element per sample; the name of the list element is the gene symbol on
which the numeric vector is a vector of expression values; the
\code{snpid} parameter is bound to a formula template, either
\code{~.} or \code{~factor(.)} expressing choice between an additive
genetic model and a general model.}
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\note{ The \code{genesym} class works in conjunction
with annotation to allow translation from HUGO symbols
to affy probe set names according to Bioconductor 
annotation environments.  The \code{getpsid} function
carries out the translation; the \code{genesym} class/wrapper
just tells us that a string will be serving the role of
a HUGO symbol and allows us to dispatch accordingly.}

\section{Notes}{Helper functions \code{SNPalleles}
and \code{rarebase} are accessors}

% ~Make other sections like Warning with \section{Warning }{....} ~
%
\seealso{
\code{\link{make_racExSet}} for construction; \code{\link[Biobase]{eSet}} ancestor
class
}
\examples{
example(make_racExSet)
dem
snpNames(dem)[1:4]
featureNames(dem)[1:4]
oneFit(dem, genesym("DDR1"), "rs4814683", lm)
ddr1 = list(DDR1=exprs(dem)["1007_s_at",])
oneFit(dem, ddr1, "rs4814683", aov)
data(chr20GGceu)
ddr1 = list(DDR1=exprs(chr20GGceu)["1007_s_at",])
summary(oneFit(chr20GGceu, ddr1, ~rs4814683+factor(rs6060535), lm))
}
\keyword{classes}
