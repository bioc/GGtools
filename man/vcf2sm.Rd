\name{vcf2sm}
\alias{vcf2sm}
\alias{show,metaVCF-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate a snp.matrix instance on the basis of a VCF (4.0) file
}
\description{
generate a snp.matrix instance on the basis of a VCF (4.0) file
}
\usage{
vcf2sm(gzpath, chrom, tabixcmd = "tabix", nmetacol = 9, verbose = FALSE, gran=10000)
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
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
an instance of \code{\link[snpMatrix]{snp.matrix-class}}
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
\dontrun{
 vref = system.file("vcf/ex.vcf.gz", package="GGtools")
 vcf2sm( vref, "20" )
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
