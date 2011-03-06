\name{imphm3_1KG_20}
\alias{imphm3_1KG_20_mA2}
\alias{imphm3_1KG_20_mA5}
\docType{data}
\title{
snpStats-generated rules from imputing from HapMap phase III loci to 1000 genomes loci -- for chromosome 20 only
}
\description{
snpStats-generated rules from imputing from HapMap phase III loci to 1000 genomes loci -- for chromosome 20 only
}
\usage{data(imphm3_1KG_20_mA2)}
\format{
  The format is:
Formal class 'snp.reg.imputation' [package "snpStats"] with 1 slots \cr
  ..@ .Data:List of 110511 \cr
  .. ..$ :List of 4 \cr
  .. .. ..$ maf         : num 0.2 \cr
  .. .. ..$ r.squared   : num 1 \cr
  .. .. ..$ snps        : chr "rs6139074" \cr
  .. .. ..$ coefficients: num [1:2] 0 1 \cr
  .. ..$ :List of 4 \cr
  .. .. ..$ maf      : num 0.117 \cr
  .. .. ..$ r.squared: num 0.892 \cr
  .. .. ..$ snps     : chr [1:3] "rs13043000" "rs17685809" "rs1935386" \cr
  .. .. ..$ hap.probs: num [1:16] 3.01e-01 6.97e-22 1.56e-02 2.36e-20 8.49e-03 ... \cr
  .. ..$ : NULL
}
\details{
Generated with snpStats 1.1.1, rules that use the ceu1kg package to
define loci and calls for 1000 genomes genotypes for CEU, to allow imputation
from the hapmap phase III loci for CEU.  The data object with suffix mA2 was generated with
setting mA=2; for suffix mA5, mA was set at 5; see \code{\link[snpStats]{snp.imputation}}
for details on this parameter, which sets the minimum number of observations required for
an LD determination to be made for SNP tagging or haplotype modeling.
}
\source{
ceuhm3 package was used to define the hapmap phase III loci; locations
derived from SNPlocs.Hsapiens.dbSNP.20090506; ceu1kg package includes
metadata and calls derived from the 1000 genomes pilot phase 1 VCF file
for CEU.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(imphm3_1KG_20_mA2)
imphm3_1KG_20_mA2[1:10]
}
\keyword{datasets}
