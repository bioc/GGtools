\name{detScreen}
\alias{detScreen}
\title{perform a screen for SNP determination of expression
using a ggExprSet instance, and produce a plot  }
\description{perform a screen for SNP determination of expression
using a ggExprSet instance, and produce a plot
}
\usage{
detScreen(gge=c20GGceu, psn="206918_s_at", chrmeta=chr20meta, chr="chr20", gran=50, gene="CPNE1", ...)
}
\arguments{
  \item{gge}{ \code{\link{ggExprSet}} instance}
  \item{psn}{ geneNames element for \code{gge}}
  \item{chrmeta}{ data.frame with chromosome location of
SNPs indexed using names of snps component of \code{gge}}
  \item{chr}{ character label for chromosome (should
be computed from \code{chrmeta} arg...)}
  \item{gran}{ deterministic selection: \code{1/gran} snps are used; if you
want to use all snps available use \code{gran}=1}
  \item{gene}{ character string with gene symbolic name for plotting purposes,
corresponds to \code{psn}}
  \item{\dots}{ additional arguments to regseq }
}
\details{
regressions for additive effect on expression of the
gene denoted by \code{psn} are run for all snps (as
selected using gran parameter) on the chromosome
}
\value{
Plots the negative log 10 p value of regression test for additive
effect of genotype on expression for the selected gene
over all snps, as well as the specific plot of regression for the
most significant snp;
Invisibly returns a list with the index of the most significant
snp, and the list of regression models for all snps (as computed
by \code{\link{regseq}}
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(c20GGceu)
data(chr20meta)
args(detScreen)
# for defaults, need to see counting message go to 2200 (granularity 50)
# rebuild compliant metadata structure
tmp = as(chr20meta,"data.frame")
chr20meta = data.frame(rsnum=rownames(tmp), tmp)
detScreen()
}
\keyword{ models }
