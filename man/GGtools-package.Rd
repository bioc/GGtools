\name{GGtools-package}
\alias{GGtools-package}
\alias{GGtools}
\alias{hg19.si.df}
\docType{package}
\title{
software and data for analyses in genetics of gene expression
}
\description{
software and data for analyses in genetics of gene expression
}
\details{
\tabular{ll}{
Package: \tab GGtools\cr
Version: \tab 4.2.26\cr
Suggests: \tab GGdata, illuminaHumanv1.db\cr
Depends: \tab R (>= 2.14), GGBase (>= 3.16.1)\cr
Imports: \tab methods, snpStats, ff, IRanges, GenomicRanges,
AnnotationDbi, Biobase, Rsamtools, bit, VariantAnnotation\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
Packaged: \tab 2012-01-18 03:39:51 UTC; stvjc\cr
Collate: \tab AllClasses.R AllGenerics.R eqtlTests.R
managers.R topFeats.R gwSnpTests.R snpsCisToGenes.R
relocate.R topSnps.R
transutils.R vcfutils.R eqtlEstimates.R alleq.R meta.R\cr
Built: \tab R 2.15.0; ; 2012-02-06 17:22:52 UTC; unix\cr
}

Index:
\preformatted{
best.cis.eQTLs          collect genewise best scoring eQTL
eqtlTests               compute association statistics between all
                        probes and SNP in an smlSet instance
eqtlTestsManager-class
                        Class '"eqtlTestsManager"'
ex                      ExpressionSet instance for illustrating
                        integrative smlSet container
getCisMap               create, using Bioconductor annotation
                        resources, a structure that enumerates SNP in
                        the vicinity of ('cis' to) genes
gwSnpTests              execute a series of tests for association
                        between genotype and expression
strMultPop              serialization of a table from Stranger's
                        multipopulation eQTL report
hg19.si.df		data frame representation of seqinfo for Homo.sapiens
			at hg19 build
}

 The package depends on GGBase, which includes additional
 infrastructure for integrative data structures and data
 filtering.

}
\author{
VJ Carey <stvjc@channing.harvard.edu>

Maintainer: VJ Carey <stvjc@channing.harvard.edu>
%~~ The author and/or maintainer of the package ~~
}
%\references{
%~~ Literature or other references for background information ~~
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
\seealso{
 \code{\link[GGBase]{getSS}} for acquiring containers for
integrative data on genetics of expression.
}
\examples{
\dontrun{
 # acquire chromosome 20 genotypes and all expression data for
 # 90 CEU samples as published at Wellcome Trust GENEVAR and
 # HapMap phase II
 c20 = getSS("GGtools", "20")
 # perform a focused eQTL search
 t1 = gwSnpTests(genesym("CPNE1")~male, c20)
 # get best hits
 topSnps(t1)
 }
}
