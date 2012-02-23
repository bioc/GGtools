\name{vcf2sm-methods}
\docType{methods}
\alias{vcf2sm-methods}
\alias{vcf2sm,TabixFile,GRanges,integer-method}
\alias{vcf2sm}
\title{  Methods for Function \code{vcf2sm} in Package \pkg{GGtools} }
\description{
A provisional transformer from tabix-indexed VCF to SnpMatrix instances
}
\section{Methods}{
\describe{

\item{\code{signature(tbxfi = "TabixFile", gr = "GRanges", nmetacol = "integer")}}{
Filters to records overlapped by \code{gr}, then
converts calls to B allele counts+1; 0 is reserved for
no-calls; \code{nmetacol} is
conventionally set to \code{9L}.

More work is needed to properly recover uncertain imputed genotypes.
See also the related \code{\link[VariantAnnotation]{readVcf}} procedure.
}
}}
\examples{
tf = TabixFile(system.file("vcf/CEU.exon.2010_09.genotypes.vcf.gz",
 package="GGtools"))
gr = GRanges("1", IRanges(1.1e6,1.15e6))
vcf2sm(tf, gr=gr, nmetacol=9L)
}
\keyword{methods}
