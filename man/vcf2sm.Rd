\name{vcf2sm}
\alias{vcf2sm}
\alias{show,metaVCF-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate a SnpMatrix instance on the basis of a VCF (4.0) file
}
\description{
generate a SnpMatrix instance on the basis of a VCF (4.0) file. NOTE: the
tabix utility must be installed and be invocable via system().
}
\usage{
vcf2sm(gzpath, chrom, tabixcmd = "tabix", nmetacol = 9, verbose = FALSE, gran=10000,
   metamax=100, makelocpref="chr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gzpath}{
string: path to a gzipped vcf file
}
  \item{chrom}{
string: chromosome for processing; use tabix -l to obtain the
list of tokens if necessary
}
  \item{tabixcmd}{
string: assumes tabix available as an executable utility; tells
the absolute path for invoking the command
}
  \item{nmetacol}{
numeric: tells number of columns used in each record as locus-level metadata
}
  \item{verbose}{
logical: if TRUE, provide processing info
}
  \item{gran}{
numeric: a report is given once every gran snp are traversed to show progress
}
  \item{metamax}{number of lines to be sniffed for metadata before real data
encountered, could be liberal}
  \item{makelocpref}{string telling what to use to construct a 
locus identifier when the id field is .; sometimes loc field is adequate and this
should be set to "".  set to "" if you see loc names with chrchr prefix.}
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
# requires tabix
chkTabix = try(system("tabix 2>&1", intern=TRUE))
if (!inherits(chkTabix, "try-error") && length(grep("Option", chkTabix))>0) {
 vref = system.file("vcf/CEU.exon.2010_09.genotypes.vcf.gz", package="GGtools")
 vcf2sm( vref, "1" )
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
