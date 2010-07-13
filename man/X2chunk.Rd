\name{X2chunk}
\alias{X2chunk}
\alias{topFeats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute numerical matrix of chisq statistics in a genomic interval;
extract features as requested
}
\description{
compute numerical matrix of chisq statistics in a genomic interval
(rows are SNP, columns are genes), or extract features
}
\usage{
X2chunk(mgr, ffind, start, end, snplocs, anno, useSym)
topFeats( probeid=NULL, sym=NULL, rsid=NULL, mgrOrCTD, ffind, anno, n=10, useSym=TRUE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mgr}{
an instance of multffManager
}
  \item{mgrOrCTD}{
an instance of multffManager or a cisTransDirector instance
}
  \item{probeid}{
character string for a probeid
}
  \item{sym}{
character string for a gene
}
  \item{rsid}{
character string for a snp id
}
  \item{ffind}{
the index of the ff structure to use (typically chromosome number)
}
  \item{start}{
left end of interval of interest
}
  \item{end}{
right end of interval of interest
}
  \item{snplocs}{
location structure for SNP (RangedData instance)
}
  \item{n}{
for topFeats, the number of features to report
}
  \item{anno}{
name of a gene annotation package resolving the identifiers
used in column names of ff matrix
}
  \item{useSym}{
logical indicating whether colnames of return should be gene symbols derived
from \code{anno}
}
}
\details{
X2chunk will obtain RAM resources for material on disk, so use with caution

Note that gene symbols may map to multiple probes.  The first hit is used
by topFeats when used with sym=.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
}


\examples{
# build an smlSet with a small set of neighboring genes
data(snpLocs20)
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
library(illuminaHumanv1.db)
gOn20 = get("20", revmap(illuminaHumanv1CHR))
gLocs = geneRanges(gOn20, "illuminaHumanv1.db")
start = 10000000
end = 13500000
g2use_inds = which(ranges(gLocs)$chr20 \%in\% IRanges(start,end))
g2use_names = gLocs[g2use_inds,]$name
h20 = hmceuB36.2021[ probeId(g2use_names), ]
h20 = h20[chrnum(20),]
sn2use_inds = which(ranges(snpLocs20)$chr20 \%in\% IRanges(start,end))
od = getwd()
setwd(tempdir())
# create the ff manager instance
dd = diagffCC(h20, gs~male)
# extract the matrix
fc = X2chunk(dd, 1, start, end, snpLocs20, "illuminaHumanv1.db")
dim(fc)
fc[1:4,1:5]
setwd(od)
heatmap(fc[1:50,], Rowv=NA, Colv=NA, scale="none")
topFeats(probeid=NULL, sym=NULL, rsid="rs6094162", mgr=dd, 1, "illuminaHumanv1.db")
topFeats(probeid=NULL, sym="MKKS", rsid=NULL, mgr=dd, 1, "illuminaHumanv1.db")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
