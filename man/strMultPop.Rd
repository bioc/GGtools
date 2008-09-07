\name{strMultPop}
\alias{strMultPop}
\docType{data}
\title{ serialization of a table from Stringer's multipopulation eQTL report }
\description{ serialization of a table from Stringer's multipopulation eQTL report }
\usage{data(strMultPop)}
\format{
  A data frame with 39649 observations on the following 12 variables.
  \describe{
    \item{\code{rsid}}{a factor with levels rs...}
    \item{\code{genesym}}{a factor with levels \code{37865} \code{39692} \code{ABC1} \code{ABCD2} \code{ABHD4} \code{ACAS2} ...}
    \item{\code{illv1pid}}{a factor with levels \code{GI_10047105-S} \code{GI_10092611-A} \code{GI_10190705-S} \code{GI_10567821-S} \code{GI_10835118-S} \code{GI_10835186-S}  ...}
    \item{\code{snpChr}}{a numeric vector}
    \item{\code{snpCoordB35}}{a numeric vector}
    \item{\code{probeMidCoorB35}}{a numeric vector}
    \item{\code{snp2probe}}{a numeric vector}
    \item{\code{minuslog10p}}{a numeric vector}
    \item{\code{adjR2}}{a numeric vector}
    \item{\code{assocGrad}}{a numeric vector}
    \item{\code{permThresh}}{a numeric vector}
    \item{\code{popSet}}{a factor with levels \code{CEU-CHB-JPT} \code{CEU-CHB-JPT-YRI} \code{CHB-JPT}}
  }
}
\details{
imported from the PDF(!) distributed by Stranger et al as supplement to 
PMID 17873874
}
\source{
PMID 17873874 supplement
}
\references{
PMID 17873874 supplement
}
\examples{
data(strMultPop)
strMultPop[1:2,]
}
\keyword{datasets}
