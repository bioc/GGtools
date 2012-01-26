\name{GGtools-package}
\alias{GGtools-package}
\alias{GGtools}
\docType{package}
\title{GGtools Package Overview}
\description{
   GGtools Package Overview
}
\details{
  This package provides facilities for analyzing relationships
  between gene expression distributions (singly or in groups)
  and SNP genotype series (chromosome-specific or genome-wide).
  The \code{\link{gwSnpTests}} method is the primary interface.

  Important data classes in use: \code{\link[GGBase]{smlSet-class}},
  \code{\link[GGBase]{gwSnpScreenResult-class}}, defined in GGBase package.

  Main data sets: \code{\link{ex}}, an ExpressionSet that can be
  linked to genotype data based on
  chromosomes 20 and 21, with genotypes for all phase II HapMap SNP
  and full expression data for 90 CEU HapMap cohort members.  To
  create this example dataset, use \code{getSS("GGtools", c("20", "21"))}

  Introductory information is available from vignettes, type
  \code{openVignette()}.

  Full listing of documented articles is available in HTML view by
  typing \code{help.start()} and selecting \code{GGtools} package from the
  \code{Packages} menu or via \code{library(help="GGtools")}.
}
\author{
   V. Carey 
}
\keyword{package}
