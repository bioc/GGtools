\name{scoresInRanges}
\alias{scoresInRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
structured survey of eqtlTestsManager to retrieve scores cis to genes
}
\description{
structured survey of eqtlTestsManager to retrieve scores cis to genes
defined through ranges
}
\usage{
scoresInRanges(mgr, geneRanges, snpRanges, applier = lapply)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mgr}{
instance of \code{\link{eqtlTestsManager-class}}
}
  \item{geneRanges}{
instance of \code{\link[GenomicRanges]{GRanges-class}} representing
gene-like regions within which SNPs will be sought for association with
expression
}
  \item{snpRanges}{
instance of \code{\link[GenomicRanges]{GRanges-class}} representing
SNP locations
}
  \item{applier}{
lapply or mclapply if relevant
}
}
%\details{
%#%%  ~~ If necessary, more details than the description above ~~
%#}
\value{
a list of lists, one element per gene range, with matrix of chisq scores
for SNP
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

\examples{
gl = GGtools:::geneLimits( anno="illuminaHumanv1.db" )
egl = GGtools:::extendGR(gl, siz=5e5)
data(hmceuB36.2021)
gl20 = names(egl[ which(seqnames(egl) == "chr20") ] )[1:20]
p2 = intersect( gl20, featureNames(hmceuB36.2021) )
curd = getwd()
setwd(tempdir())
if (file.exists("foo"))system("rm -rf foo")
ee = eqtlTests( hmceuB36.2021[ probeId(p2), ], ~1 )
library(ceu1kg)
data(ceu1kgMeta_20)
ceu1kgMeta_20 = updateObject(ceu1kgMeta_20)
cc = scoresInRanges( ee, egl[1:10], ceu1kgMeta_20 )
sapply(cc, sapply, length)
sapply(cc, sapply, max)
setwd(curd)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
