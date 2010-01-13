\name{masterSnps}
\alias{masterSnps}
\title{ visualize a multiGwSnpScreenResult}
\description{ visualize a multiGwSnpScreenResult}
\usage{
masterSnps(mgw, n = 50, auto = TRUE, orgdb = "org.Hs.eg.db", minl10 = 5, 
gstart = 0, gend = 3e+09,
genomesize = 3e+09, pcex = 1, pal = rainbow(20), numxax=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mgw}{ a multiGwSnpScreenResult, for example from
gwSnpTests with a GeneSet on lhs of formula }
  \item{n}{ number of best snps to retain per gene}
  \item{auto}{ restrict attention to autosomes? }
  \item{orgdb}{ an annotation library like org.Hs.eg.db }
  \item{minl10}{ threshold of -log10 p above which we keep SNP for plotting}
  \item{gstart}{position at which genome-wide SNP locations begin}
  \item{gend}{position at which genome-wide SNP locations end}
  \item{genomesize}{ number of bases over which plotting will be conducted
(e.g., ylim=c(0, genomesize) ) }
  \item{pcex}{ cex setting for pch of plot }
  \item{pal}{ a palette to differentiate gene coloring }
  \item{numxax}{ logical: if TRUE, x axis labels genomic coordinates,
otherwise chromosome }
  \item{\dots}{ args passed to plot() }
}
\details{
experimental display with snp location as ordinate and
gene location as abscissa -- point plotted if snp is
associated with gene at p smaller than the threshold specified
}
\value{
a list with self-describing elements
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\examples{
if (require("GGdata")) {
 data(fheadFits)
 mm7 = masterSnps(fheadFits, minl10=7, pal=rainbow(10))
 }
}
\keyword{ models }
