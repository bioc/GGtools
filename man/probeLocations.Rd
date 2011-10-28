\name{probeLocations}
\alias{probeLocations}
\alias{probeChromosomes}
\alias{probeSequences}
\alias{snpLocations}
\alias{proximityList}
\alias{restrictProbesToChrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
utilities for annotation acquisition with smlSet instances
}
\description{
utilities for annotation acquisition with smlSet instances
}
\usage{
probeLocations(sms, extend=0)
probeSequences(sms)
probeChromosomes(sms)
snpLocations(sms, snpLocGRanges, grsnpid = "RefSNP_id")
proximityList(sms, smlind=1, snpLocGRanges, grsnpid = "RefSNP_id", probeLocExtend=0,
  glocTransform=function(x)x)
restrictProbesToChrom(smlSet, chrom)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sms}{
instance of \code{\link[GGBase]{smlSet-class}}, which should have
an annotation slot correctly identifying an AnnotationDbi-compliant
annotation package for expression probes
}
  \item{smlSet}{see \code{sms} above}
  \item{extend}{numeric value to extend regions by, using algebra of + for IRanges instances}
  \item{snpLocGRanges}{an instance of \code{\link[GenomicRanges]{GRanges-class}} with numeric tokens for dbSNP identifiers as used in SNPlocs.* packages.
The identifiers are in elementMetadata(snpLocGRanges)[[grsnpid]]
}
  \item{grsnpid}{a character string naming the elementMetadata element
holding the SNP identifiers, as returned by a getSNPlocs with
as.GRanges=TRUE}
  \item{smlind}{\code{proximityList} will work for only one
element of \code{smList(sms)}; this argument picks it out.}
  \item{probeLocExtend}{numeric, by which we extend (using ranges()+extend)
the ranges of the resulting
GRanges for probe locations}
  \item{glocTransform}{alters the \code{probeLocations(sms)} GRanges
to define SNP proximity for inclusion in gene-snp association list}
  \item{chrom}{chromosome token}
} % end args
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{probeLocations} and \code{snpLocations} return
GRanges instances.  names() of the probeLocations result
gives the probe identifiers.  SNP identifiers are in
elementMetadata.

probeSequences() works when a PROBESEQUENCE custom mapping is available
in the AnnotationDbi compliant package identified in \code{annotation(sms)}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(smlSet.example)
library(illuminaHumanv1.db)
probeLocations(smlSet.example)
spac = "SNPlocs.Hsapiens.dbSNP.20100427"
if (spac \%in\% (IP <- installed.packages()[,1])) {
 library(spac, character.only=TRUE)
 c20 = getSNPlocs("ch20", as.GRanges=TRUE)
 seqlevels(c20)[20] = "20"
 if (!exists("hmceuB36.2021")) hmceuB36.2021 <- getSS("GGtools", c("20", "21"))
 h20 = hmceuB36.2021[ chrnum("20"), ] # restrict smList
 pcc = probeChromosomes(h20)
 p20 = featureNames(h20)[ which(pcc == "20") ] # for NA
 h20 = h20[ probeId(p20[1:20]), ]
 h20 = dropMonomorphies(h20)
 pl = proximityList( h20, 1, c20)
 cat("count SNPs in gene region\n")
 print(sapply(pl, length))
 pl2 = proximityList( h20, 1, c20, glocTransform = function(x) {
    ranges(x) = ranges(x)+50000 # extend gene interval 50kb on each end
    x})
 pl3 = proximityList( h20, 1, c20, probeLocExtend=50000)
 cat("count SNPs in gene region extended by 50k at each end\n")
 print( sapply(pl2, length))
 td = tempdir()
 curd = getwd()
 setwd(td)
 applier = lapply
 if ("multicore" \%in\% IP) {
    library(multicore)
    applier=mclapply
    }
 et1 = eqtlTests( h20, ~male, geneApply=applier )
 sco1 = lapply(1:length(pl), function(x) et1[ rsid(pl[[x]]), probeId(names(pl)[x]) ])
 print(sapply(sco1, sapply, max))
 sco2 = lapply(1:length(pl2), function(x) et1[ rsid(pl2[[x]]), probeId(names(pl2)[x]) ])
 print(sapply(sco2, sapply, max))
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
