\name{m20}
\alias{m20}
\docType{data}
\title{
snpMatrix2 (1.1.1) with imputed genotypes for 110
HapMap phase III samples from CEU population
}
\description{
snpMatrix2 (1.1.1) with imputed genotypes for 110
HapMap phase III samples from CEU population
}
\usage{data(m20)}
\format{
  The format is:
Formal class 'snp.matrix' [package "snpMatrix2"] with 1 slots \cr
  ..@ .Data: raw [1:110, 1:190473] 03 03 03 03 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:110] "NA06984" "NA06989" "NA12340" "NA12341" ... \cr
  .. .. ..$ : chr [1:190473] "rs6078030" "rs4814683" "rs34147676" "rs6139074" ... \cr
}
\details{
results of MACH applied by Blanca Himes of Channing Laboratory,
leading to an mlprob file read with read.mach()
}
\source{
The HapMap phase III genotypes were obtained as
\code{hapmap3_r2_b36_fwd.CEU.qc.poly.[ped/map]} as distributed at
hapmap.org
}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(m20)
}
\keyword{datasets}
