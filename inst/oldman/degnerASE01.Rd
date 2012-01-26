\name{degnerASE01}
\alias{degnerASE01}
\docType{data}
\title{
transcription of a table from a paper by Degner et al
}
\description{
transcription of a table from a paper by Degner et al, involving identification of
genes with allele-specific expression discovered by RNA-seq
}
\usage{data(degnerASE01)}
\format{
  A data frame with 55 observations on the following 10 variables.
  \describe{
    \item{\code{rsnum}}{a factor with levels \code{rs10266655} \code{rs1042448} \code{rs1046747} \code{rs1047469} \code{rs1059307} \code{rs1060915} \code{rs11009147} \code{rs1127326} \code{rs11376} \code{rs11570126} \code{rs11578} \code{rs1158} \code{rs13306758} \code{rs13309} \code{rs16952692} \code{rs17014852} \code{rs17459} \code{rs1879182} \code{rs2070924} \code{rs2071888} \code{rs2089910} \code{rs2234978} \code{rs2271920} \code{rs2530680} \code{rs3025040} \code{rs3170545} \code{rs325400} \code{rs368116} \code{rs3819946} \code{rs3871984} \code{rs4784800} \code{rs4982685} \code{rs558018} \code{rs6568} \code{rs6682136} \code{rs6890805} \code{rs7046} \code{rs705} \code{rs7121} \code{rs7141712} \code{rs7192} \code{rs7695} \code{rs7739387} \code{rs8023358} \code{rs8084} \code{rs8429} \code{rs8647} \code{rs8905} \code{rs9038} \code{rs916974}}
    \item{\code{refreads}}{a numeric vector}
    \item{\code{nonrefreads}}{a numeric vector}
    \item{\code{miscall}}{a numeric vector}
    \item{\code{chr}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr22} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9}}
    \item{\code{loc}}{a numeric vector}
    \item{\code{gene}}{a factor with levels \code{ADAR} \code{ADPGK} \code{AKAP2} \code{AP4M1} \code{ATF5} \code{BIN1} \code{BRCA1} \code{C6orf106} \code{CCL22} \code{CD59} \code{CRYZ} \code{DFNA5} \code{ENSA} \code{FAS} \code{GNAS} \code{GYPC} \code{HLA-DPB1} \code{HLA-DRA} \code{HMMR} \code{ITGB1} \code{LSP1} \code{MADD} \code{MARK3} \code{ME2} \code{MEF2A} \code{MGAT1} \code{MRPL52} \code{MTMR2} \code{NF2} \code{NIN} \code{NUP62} \code{OAS2} \code{PALM2-AKAP2} \code{PIP4K2A} \code{PRKAR1A} \code{PTK2B} \code{SAR1A} \code{SEC22B} \code{SEMA4A} \code{SEPT9} \code{SLC2A1} \code{SNHG5} \code{SNURF/SNRPN} \code{STX16} \code{TAF6} \code{TAPBP} \code{VEGFA}}
    \item{\code{indiv}}{a factor with levels \code{GM19238} \code{GM19239}}
    \item{\code{eqtl}}{a factor with levels \code{} \code{Yes}}
    \item{\code{imprint}}{a logical vector}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Effect of read-mapping biases on detecting allele-specific expression from RNA-sequencing data.
Jacob F. Degner 1,3,, John C. Marioni 1,, Athma A. Pai 1, Joseph K. Pickrell 1, Everlyne Nkadori 1,2, Yoav Gilad 1, and Jonathan K. Pritchard 1,2, Bioinformatics 2009.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(degnerASE01)
degnerASE01[1:4,]
## maybe str(degnerASE01) ; plot(degnerASE01) ...
}
\keyword{datasets}
