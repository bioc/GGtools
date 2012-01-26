\name{best.cis.eQTLs}
\alias{best.cis.eQTLs}
\alias{mcwBestCis-class}
\alias{chromsUsed}
\alias{chromsUsed,mcwBestCis-method}
\alias{fullreport}
\alias{fullreport,mcwBestCis,missing-method}
\alias{fullreport,mcwBestCis,character-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
collect genewise best scoring eQTL
}
\description{
collect genewise best scoring eQTL
}
\usage{
best.cis.eQTLs(smpack = "GGdata", rhs = ~1, 
  folderstem = "cisScratch", radius = 50000, 
  chrnames = as.character(1:22), 
  smchrpref = "", gchrpref = "", schrpref = "ch", 
  geneApply = lapply, geneannopk = "illuminaHumanv1.db", 
  snpannopk = "SNPlocs.Hsapiens.dbSNP.20100427", 
  smFilter = function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97), nperm = 2, cisMapList = NULL)

chromsUsed(x)

fullreport(x, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smpack}{
character string naming a package to which \code{\link[GGBase]{getSS}}
can be applied to extract \code{\link[GGBase]{smlSet-class}} instances
}
  \item{rhs}{
R model formula, with no dependent variable, that will be used with
\code{\link[snpStats]{snp.rhs.tests}} to adjust GWAS tests for
each expression probe
}
  \item{folderstem}{
prefix of the folder name to be used to hold ff archives of
test results
}
  \item{radius}{
coding extent of each gene will be extended in both directions by
\code{radius} bases, and only SNP within these limits are
used for selecting best hits for the gene
}
  \item{chrnames}{
character vector of chromosome identifiers, to be manipulated
for certain query resolutions by the following parameters
}
  \item{smchrpref}{
prefix to convert \code{chrnames} into appropriate tokens for
indexing \code{smlSet} elements as collected from the
package named by parameter \code{smpack}
}
  \item{gchrpref}{
prefix to convert \code{chrnames} into appropriate tokens for
obtaining gene metadata; in future this may need to be a string
transformation function
}
  \item{schrpref}{
prefix to convert \code{chrnames} into appropriate tokens for
use with \code{getSNPlocs} for the SNP location information
package identified in \code{snpannopack} parameter below
}
  \item{geneApply}{
an lapply like function, defaults to \code{lapply}
}
  \item{geneannopk}{
character string, name of annotation package that annotates probe identifiers
}
  \item{snpannopk}{
character string, name of SNPlocs.Hsapiens.dbSNP.* package for
obtaining 
}
  \item{smFilter}{
function accepting and returning an \code{\link[GGBase]{smlSet-class}}
instance
}
  \item{nperm}{
number of permutations to be used for plug-in FDR computation
}
  \item{cisMapList}{
if supplied, must be a list with one element per element of
\code{chrnames}; each element is a list with one element per
probe resident on the corresponding chromosome, and the
probe element's content is a vector of SNP identifiers that are
deemed `cis'.  This list will be computed using annotation resources
if this argument is left missing.
}
  \item{x}{instance of \code{mcwBestCis}}
  \item{type}{character, either 'data.frame' or 'GRanges'}
}
\details{
\code{geneApply} can be set to \code{parallel::mclapply}, for example,
in a multicore context.

\code{mcwBestCis} stands for 'multi-chromosome-wide best cis'
eQTL report container.
}
\value{
an instance of \code{\linkS4class{mcwBestCis}}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
getClass("mcwBestCis")
\dontrun{
best.cis.eQTLs(chrnames="20")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
