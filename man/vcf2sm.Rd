\name{vcf2sm}
\alias{vcf2sm}
\alias{show,metaVCF-method}
\alias{vcf2sm,TabixFile,GRanges,integer-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate a SnpMatrix instance on the basis of a VCF (4.0) file
}
\description{
generate a SnpMatrix instance on the basis of a VCF (4.0) file. 
}
\usage{
vcf2sm(tbxfi, ..., gr, nmetacol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tbxfi}{instance of \code{\link[Rsamtools]{TabixFile-class}}}
  \item{\dots}{not used}
  \item{gr}{instance of \code{\link[GenomicRanges]{GRanges-class}}}
  \item{nmetacol}{numeric: tells number of columns used in each record as locus-level metadata }
}
\details{
This function is relevant only for diallelic SNP.  If any base call is denoted `.', the associated genotype is
set to missing (raw 0), even if the nonmissing call is ALT, implying at least one ALT.  
}

\value{
an instance of \code{\link[snpStats]{SnpMatrix-class}}
}
\references{
\url{http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf4.0}
}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# SRC: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/exon/CEU.exon.2010_03.genotypes.vcf.gz
 vref = system.file("vcf/CEU.exon.2010_09.genotypes.vcf.gz", package="GGtools")
 gg = GenomicRanges::GRanges(seqnames="1", IRanges::IRanges(10e6,20e6))
 vcf2sm(Rsamtools::TabixFile(vref), gr=gg, nmetacol=9L)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
