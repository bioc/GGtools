\name{snpScreen}
\alias{snpScreenResult-class}
\alias{twSnpScreenResult-class}
\alias{show,snpScreenResult-method}
\alias{snpScreen,racExSet,snpMeta,genesym,formula,function,numeric-method}
\alias{snpScreen,racExSet,snpMeta,genesym,formula,function,missing-method}
\alias{snpScreen,racExSet,snpMeta,genesym,formula,GGfitter-method}
\alias{twSnpScreen,racExSet,snpMeta,formula,GGfitter-method}
\alias{twSnpScreen}
\alias{snpScreen}
\alias{extract_p}
\alias{plot_mlp}
\title{compute model fits over a sequence of SNPs  }
\description{
compute model fits over a sequence of SNPs
}
\usage{
snpScreen(racExSet, snpMeta, gene, formTemplate, fitter, gran)
twSnpScreen(racExSet, snpMeta, formTemplate, fitter)
extract_p(ssr)
plot_mlp(ssr, snpMeta, gchr=NULL, geneLocDF=NULL, ps = NULL, pch = 20, cex = 0.5, local=FALSE, 
  plotf=smoothScatter, organism="human")
}
\arguments{
  \item{racExSet}{ instance of \code{\link{racExSet-class}}}
  \item{snpMeta}{ instance of \code{\link{snpMeta-class}}}
  \item{gene}{ instance of \code{\link{genesym-class}} class identifying
the expression phenotype to be regarded as dependent variable}
  \item{gchr}{ string stating chromosome on which gene is resident; looked up in geneLocDF if not supplied}
  \item{geneLocDF}{ data frame that has one row for each gene symbol, in column `gene', and
     the chromosome on which it lives, in column `chr'}
  \item{formTemplate}{ a formula having form \code{~.} or
\code{~factor(.)}, literally, to specify additive or
nonadditive models for effects of rare allele copy number}
  \item{fitter}{ R fitting function that can work with formulas
and data frames, for example \code{\link{lm}}, or, if not a
fitting function, a special object in the package, either
\code{fastHETfitter}, or \code{fastAGMfitter}, see details below}
  \item{gran}{ a numeric convenience parameter; if not 1, SNPs will
be deterministically selected with frequency 1/\code{gran} along the
linear sequence in the racExSet racAssay structure }
  \item{\dots}{ \dots -- not in use}
  \item{ssr}{\code{snpScreenResult} instance}
  \item{ps}{a string with a probe set identifier}
  \item{pch}{for passage to \code{\link[geneplotter]{scatterSmooth}}}
  \item{cex}{for passage to \code{\link[geneplotter]{scatterSmooth}}}
  \item{local}{by default (local==FALSE) the plot is made over the
entire chromosome; if local==TRUE, the plot is made over the
segment of the chromosome within which the snps screened lie;
this only makes a difference if the racExSet used in the screen
has been SNP-filtered relative to the HapMap SNP set. }
  \item{plotf}{function to use for rendering -- with really
high-density genotyping, the default works well; in sparse cases,
use \code{\link{plot}}.}
  \item{organism}{string used for plot annotation}
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

For twSnpScreen (transcriptome-wide): all genes in the expression array are
checked against all SNPs.  An instance of twSnpScreenResult is
returned.

\code{fastAGMfitter} is an object that wraps a
call to a C routine for simple least square fitting
of an additive genetic model.  \code{fastHETfitter} will call
C code that compares heterozygous
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
dem = exclMono(dem)
snpNames(dem)[1:4]
featureNames(dem)[1:4]
data(chr20GGdem)
data(chr20meta)
data(geneLocs_hsa)
scr1 = snpScreen(dem, chr20meta, genesym("DDR1"), ~., lm,  gran=100 )
scr1[[1]]
scr2 = snpScreen(dem, chr20meta, genesym("DDR1"), ~factor(.), lm,  gran=200 )
scr2[[1]]
plot_mlp(scr1, chr20meta, geneLocDF=geneLocs_hsa)
chr20GGdem = exclMono(chr20GGdem)
ut = unix.time(scr2 <- snpScreen(chr20GGdem, chr20meta, genesym("CPNE1"), ~., fastAGMfitter))
ut
scr2
plot_mlp(scr2, chr20meta, geneLocDF=geneLocs_hsa)
#
# here we work on a WebQTL computation
#
# get the expr+genotype data
data(gse2031GG)
# get a provisional snp metadata structure
data(INB34snpMeta)
# run a screen for Erdr1
ss = snpScreen(gse2031GG, INB34snpMeta, genesym("Erdr1"), ~., fastAGMfitter)
plot_mlp(ss, INB34snpMeta, gchr="all", plotf=plot, organism="mouse" )
#
# following will throw warnings, so dontrun but paste by hand
#
\dontrun{
 someg = c("DDR1", "RFC2", "HSPA6", "PAX8", "GUCA1A")
 somesnps = c("rs6088779", "rs6120881", "rs7346331", "rs6088780", "rs2425024", 
 "rs6088781", "rs6141540", "rs11167262", "rs6119594", "rs1555322", 
 "rs2425027", "rs6142341", "rs6058215", "rs6060331", "rs6088785", 
 "rs6087694", "rs12106107", "rs2425032", "rs2425037", "rs2275274")
 
 somesnps = somesnps[c(5L, 8L, 10L, 11L, 14L, 18L, 19L, 20L)]
 trim = chr20GGdem[ genesym(someg[1:3]), ]
 trim = trim[ snpID(somesnps), ]
 lktw = twSnpScreen( trim, chr20meta, ~., fastAGMfitter)
 names(lktw)
 lktw[[1]]
}

}
\keyword{ models }
