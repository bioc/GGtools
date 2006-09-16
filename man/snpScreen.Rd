\name{snpScreen}
\alias{snpScreenResult-class}
\alias{show,snpScreenResult-method}
\alias{snpScreen,racExSet,snpMeta,genesym,formula,function,numeric-method}
\alias{snpScreen}
\alias{extract_p}
\alias{fastAGM}
\alias{fastHET}
\alias{plot_mlp}
\title{compute model fits over a sequence of SNPs  }
\description{
compute model fits over a sequence of SNPs
}
\usage{
snpScreen(racExSet, snpMeta, gene, formTemplate, fitter, gran, ...)
extract_p(ssr)
plot_mlp(ssr, snpMeta, ps = NULL, pch = 20, cex = 0.5, local=FALSE)
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
and data frames, for example \code{\link{lm}}, or, if not a
fitting function, a special function in the package, either
\code{fastHET}, or \code{fastAGM}, see details below}
  \item{gran}{ a numeric convenience parameter; if not 1, SNPs will
be deterministically selected with frequency 1/\code{gran} along the
linear sequence in the racExSet racAssay structure }
  \item{\dots}{ \dots -- not in use}
  \item{ssr}{\code{snpScreenResult} instance}
  \item{ps}{a string with a probe set identifier}
  \item{pch}{for passage to \code{\link[geneplotter]{scatterSmooth}}}
  \item{cex}{for passage to \code{\link[geneplotter]{scatterSmooth}}}
  \item{local}{by default (local=FALSE) the plot is made over the
entire chromosome; if local==TRUE, the plot is made over the
segment of the chromosome within which the snps screened lie;
this only makes a difference if the racExSet used in the screen
has been SNP-filtered relative to the HapMap SNP set. }
}
\details{
for snpScreen:
If \code{options()$verbose == TRUE} then every 100th index in the
vector of snps is printed to stdout to show rate of progress.

\code{snpScreenResult} is a container for relevant information
about a screen, including a list of fit objects.

For result processing, many SNPs have no variation in observed
samples and statistical tests of association are indeterminate.
NAs will be returned for tests on such SNPs.

\code{fastAGM} is a C routine for simple least square fitting
of an additive genetic model.  \code{fastHET} will compare heterozygous
to homozygous.

\code{plot_mlp} returns a list of x and y values for the
plotted points.
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
example(make_racExSet)
dem
snpNames(dem)[1:4]
featureNames(dem)[1:4]
data(chr20meta)
scr1 = snpScreen(dem, chr20meta, genesym("DDR1"), ~., lm,  gran=100 )
scr1[[1]]
scr2 = snpScreen(dem, chr20meta, genesym("DDR1"), ~factor(.), lm,  gran=200 )
scr2[[1]]
plot_mlp(scr1, chr20meta)
chr20GGdem = exclMono(chr20GGdem)
ut = unix.time(scr2 <- snpScreen(chr20GGdem, chr20meta, genesym("CPNE1"), ~factor(.), fastAGM, 1))
ut
scr2
plot_mlp(scr2, chr20meta)
}
\keyword{ models }
