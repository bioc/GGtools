\name{gwSnpScreen-methods}
\docType{methods}
\alias{gwSnpScreen-methods}
\alias{gwSnpScreen}
\alias{cwSnpScreen}
\alias{cnumOrMissing-class}
\alias{gwSnpScreen,genesym,smlSet-method}
\alias{gwSnpScreen,genesym,smlSet,cnumOrMissing-method}
\alias{gwSnpScreen,GeneSet,smlSet,cnumOrMissing-method}
\alias{gwSnpScreen,probeId,smlSet,cnumOrMissing-method}

\title{genome-wide test for association between SNP and phenotype}
\description{genome-wide test for association between SNP and phenotype}
\section{Methods}{
\describe{
\item{sym = "genesym", sms = "smlSet", cnum="cnumOrMissing"}{ \code{sym} is an instance
of \code{genesym} class (typically just a wrapped character string
giving a SYMBOL value for the expression annotation package in use, or
a GeneSet instance, see below);
\code{sms} is an integrative instance of \code{smlSet} containing
snp.matrix genotype data in a list and coordinated expression data.
Restriction of analysis to specific chromosomes can be accomplished by
passing an instance of \code{chrnum} as the third argument.

The \code{gwSnpScreen} computes the snp.rhs.tests from snpMatrix
for all snp against the specified expression phenotype.

If \code{sym} is an instance of the GSEABase GeneSet class,
a list of gwSnpScreenResults is packaged in a multiGwSnpScreenResult
instance.

A focused selection of a probe, for those genes mapping to multiple
probes, is facilitated by casting the sym parameter input to
class probeId.  The string cast in this way is used to index into
the rows of exprs(sms).
}}}
\keyword{methods}
\examples{
data(hmceuB36.2021)
objects("package:GGBase")
g1 = gwSnpScreen( genesym("CPNE1"), hmceuB36.2021, chrnum(20) )
g1
plot(g1)
}
\keyword{ models }
