\name{snpScreen}
\alias{snpScreenResult-class}
\alias{show,snpScreenResult-method}
\alias{snpScreen,racExSet,snpMeta,genesym,formula,function,numeric-method}
\alias{snpScreen}
\title{compute model fits over a sequence of SNPs  }
\description{
compute model fits over a sequence of SNPs
}
\usage{
snpScreen(racExSet, snpMeta, gene, formTemplate, fitter, gran, ...)
}
\arguments{
  \item{racExSet}{ instance of \code{\link{racExSet-class}}}
  \item{snpMeta}{ instance of \code{\link[GGdata]{snpMeta}}}
  \item{gene}{ instance of \code{\link{genesym}} class identifying
the expression phenotype to be regarded as dependent variable}
  \item{formTemplate}{ a formula having form \code{~.} or
\code{~factor(.)}, literally, to specify additive or
nonadditive models for effects of rare allele copy number}
  \item{fitter}{ R fitting function that can work with formulas
and data frames, for example \code{\link{lm}}}
  \item{gran}{ a numeric convenience parameter; if not 1, SNPs will
be deterministically selected with frequency 1/\code{gran} along the
linear sequence in the racExSet racAssay structure }
  \item{\dots}{ \dots -- not in use}
}
\details{
If \code{options()$verbose == TRUE} then every 100th index in the
vector of snps is printed to stdout to show rate of progress.

\code{snpScreenResult} is a container for relevant information
about a screen, including a list of fit objects.
}
\value{
creates a list of model fit results; \code{\link{try}} is used to
allow failure of fit (e.g., \code{\link{lm}} may fail if a singular
model matrix is computed
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
# as in make_racExSet
data(litsn)
data(litex)
data(pdd)
dem = make_racExSet( exprs=litex, racs=t(litsn), pd=pdd, mi=new("MIAME"),
  anno="hgfocus")
dem
snpNames(dem)[1:4]
featureNames(dem)[1:4]
scr1 = snpScreen(dem, chr20meta, genesym("DDR1"), ~., lm, gran=100 )
scr2 = snpScreen(dem, chr20meta, genesym("DDR1"), ~factor(.), lm,  gran=100 )
scr1[[1]]
scr2[[1]]
}
\keyword{ models }
