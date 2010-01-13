\name{hbTests-methods}
\docType{methods}
\alias{hbTests-methods}
\alias{hbTests}
\alias{hbTests,genesym,smlSet,chrnum,numeric,numeric-method}
\title{ haplotype-block based tests for structured expression variation }
\description{ haplotype-block based tests for structured expression variation }
\section{Methods}{
\describe{

\item{fmla = "genesym", sms = "smlSet", cnum = "chrnum", rsid = "numeric", rad = "numeric"}{ 
    expression data for gene identified by \code{genesym} is extracted from \code{sms},
    and genotype data within \code{rad} base pairs of \code{rsid} are obtained and
    processed by \code{\link[mapLD]{mapLD}} to define haplotype blocks and the SNP tagging
    these blocks.  Score tests are then computed for the association of expression
    of the gene identified by \code{genesym} with haplotype copy number (additive model by
    default, but options captured by ... are passed to \code{\link[haplo.stats]{haplo.score}}.)
  }
}}
\examples{
library(GGtools)
data(hmceuB36.2021)
hmFou = hmceuB36.2021[, hmceuB36.2021$isFounder==TRUE]
hh = hbTests(genesym("CPNE1"), hmFou, chrnum(20), 33600000, 2e4 )
hh
pvals(hh)
plot(locs(hh), -log10(pvals(hh)))
hscores(hh)[[which.min(pvals(hh))]]
}
\keyword{methods}
