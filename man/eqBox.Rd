\name{eqBox}
\alias{eqBox}
\alias{eqDesc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
descriptive plot of expression against genotype for cisAssoc
results
}
\description{
descriptive plot of expression against genotype for cisAssoc
results
}
\usage{
eqBox(gene, snp, se, tf, radius=1e6, genome="hg19", ...)
eqDesc(gene, snp, se, tf, radius=1e6, genome="hg19", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene}{
identifier of gene in SummarizedExperiment se, must be present
in rownames(se)
}
  \item{snp}{
identifier of variant in VCF referenced by tf
}
  \item{se}{
SummarizedExperiment instance
}
  \item{tf}{
TabixFile reference for a VCF file that has been bgzipped
and tabix-indexed
}
 \item{radius}{
in order to limit the VCF import, we filter variants to those within a radius
around the selected gene -- it is assumed that the selected snp will
exist in that region (we can't extract SNP by name from vcf...)
}
\item{genome}{a simple genome identifier tag}
  \item{\dots}{
(not used with eqDesc)
will embellish plot; xlab and ylab already taken care of
}
}
\details{
The chromosome names in the VCF and the seqlevelsStyle of the Summarized
Experiment must match.  
}
\value{
eqBox produces a boxplot for all categories (including NA) of genotype.
eqDesc tabulates the genotype calls in categories.
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (require(VariantAnnotation)) {
 data(lgeu)
 lgeue = lgeu[,which(lgeu$popcode=="CEU")]
 tf20 = TabixFile(system.file("vcf/c20exch.vcf.gz", package="GGtools"))
 eqBox( "ENSG00000126005.10", "rs2425038", lgeue, tf20 )
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
