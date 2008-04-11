\name{gwSnpScreen-methods}
\docType{methods}
\alias{gwSnpScreen-methods}
\alias{gwSnpScreen}
\alias{cwSnpScreen}
\alias{cnumOrMissing-class}
\alias{gwSnpScreen,genesym,smlSet-method}
\alias{gwSnpScreen,genesym,smlSet,cnumOrMissing-method}
\title{genome-wide test for association between SNP and phenotype}
\description{genome-wide test for association between SNP and phenotype}
\section{Methods}{
\describe{
\item{sym = "genesym", sms = "smlSet", cnum="cnumOrMissing"}{ \code{sym} is an instance
of \code{genesym} class (typically just a wrapped character string
giving a SYMBOL value for the expression annotation package in use);
\code{sms} is an integrative instance of \code{smlSet} containing
snp.matrix genotype data in a list and coordinated expression data.
Restriction of analysis to specific chromosomes can be accomplished by
passing an instance of \code{chrnum} as the third argument.

The \code{gwSnpScreen} computes the snp.rhs.tests from snpMatrix
for all snp against the specified expression phenotype.
}}}
\keyword{methods}
\examples{
data(hmceuB36.2021)
library(GGdata)
objects("package:GGdata")
g1 = gwSnpScreen( genesym("CPNE1"), hmceuB36.2021, chrnum(20) )
g1
plot(g1)
}
\keyword{ models }
