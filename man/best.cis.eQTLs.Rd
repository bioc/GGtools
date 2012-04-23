\name{best.cis.eQTLs}
\alias{best.cis.eQTLs}
\alias{meta.best.cis.eQTLs}
\alias{meta.All.cis.eQTLs}
\alias{All.cis.eQTLs}
\alias{mcwBestCis-class}
\alias{allSigCis-class}
\alias{show,allSigCis-method}
\alias{show,mcwBestCis-method}
\alias{show,cwBestCis-method}
\alias{show,allCigCis-method}
\alias{chromsUsed}
\alias{fdr}
\alias{chromsUsed,mcwBestCis-method}
\alias{fullreport}
\alias{fullreport,mcwBestCis,missing-method}
\alias{fullreport,mcwBestCis,character-method}
\alias{getAll}
\alias{getBest}
\alias{getCall}

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
  shortfac = 100,
  chrnames = as.character(1:22), 
  smchrpref = "", gchrpref = "", schrpref = "ch", 
  geneApply = lapply, geneannopk = "illuminaHumanv1.db", 
  snpannopk = "SNPlocs.Hsapiens.dbSNP.20100427", 
  smFilter = function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97), nperm = 2,
  useME=FALSE, excludeRadius=NULL)

All.cis.eQTLs(maxfdr = 0.05, inbestcis = NULL, smpack = "GGdata", 
    rhs = ~1, folderstem = "cisScratch", radius = 50000, 
    shortfac = 100,
    chrnames = as.character(1:22), 
    smchrpref = "", gchrpref = "", schrpref = "ch", 
    geneApply = lapply, geneannopk = "illuminaHumanv1.db", 
    snpannopk = "SNPlocs.Hsapiens.dbSNP.20100427", 
    smFilter4cis = function(x) nsFilter(MAFfilter(clipPCs(x, 
        1:10), lower = 0.05), var.cutoff = 0.85), 
    smFilter4all = function(x) MAFfilter(clipPCs(x, 
        1:10), lower = 0.05), 
    nperm = 2)

meta.best.cis.eQTLs(smpackvec = c("GGdata", "hmyriB36"), rhslist = list(~1, 
    ~1), folderstem = "cisScratch", radius = 50000, shortfac = 100, 
    chrnames = as.character(1:22), smchrpref = "", gchrpref = "", 
    schrpref = "ch", geneApply = lapply, geneannopk = "illuminaHumanv1.db", 
    snpannopk = "SNPlocs.Hsapiens.dbSNP.20100427", smFilterList = list( 
   function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97), 
   function(x) nsFilter(MAFfilter(x, lower = 0.05), var.cutoff = 0.97) ),
 nperm = 2) 

meta.All.cis.eQTLs(maxfdr = 0.05, inbestcis = NULL, smpackvec = c("GGdata", "hmyriB36"),
  rhslist = list(~1, ~1), folderstem = "cisScratch",
  radius = 50000, shortfac=100, chrnames = as.character(1:22), smchrpref = "",
  gchrpref = "", schrpref = "ch", geneApply = lapply,
  geneannopk = "illuminaHumanv1.db",
  snpannopk = "SNPlocs.Hsapiens.dbSNP.20100427",
  smFilter4cis = function(x) nsFilter(MAFfilter(clipPCs(x, 1:10),
  lower = 0.05), var.cutoff = 0.85),
  smFilter4all = function(x) MAFfilter(clipPCs(x, 1:10),
  lower = 0.05),
  nperm = 2)


chromsUsed(x)

fdr(x)

fullreport(x, type)

getAll(x)

getBest(x)

getCall(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smpack}{
character string naming a package to which \code{\link[GGBase]{getSS}}
can be applied to extract \code{\link[GGBase]{smlSet-class}} instances
}
  \item{smpackvec}{vector of character strings naming packages
that can be used as \code{smpack} values in a series of
\code{best.cis.eQTLs} calls, one per population for meta-analysis}
  \item{rhs}{
R model formula, with no dependent variable, that will be used with
\code{\link[snpStats]{snp.rhs.tests}} to adjust GWAS tests for
each expression probe
}
  \item{rhslist}{a list of model formulae to be used as \code{rhs}
in a series of \code{best.cis.eQTLs} calls, one per population for meta-analysis}
  \item{folderstem}{
prefix of the folder name to be used to hold ff archives of
test results
}
  \item{radius}{
coding extent of each gene will be extended in both directions by
\code{radius} bases, and only SNP within these limits are
used for selecting best hits for the gene
}
  \item{shortfac}{a numeric that will scale up the
chi-squared statistic before it is converted to 
short integer for storage in ff array}
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
  \item{smFilterList}{
list of functions, one element per smlSet package used in meta analysis,
 accepting and returning an \code{\link[GGBase]{smlSet-class}}
instance
}
  \item{nperm}{
number of permutations to be used for plug-in FDR computation}
  \item{useME}{logical; if TRUE, use the rudimentary interface
to the MatrixEQTL package from A. Shabalin on CRAN }
\item{maxfdr}{Used in \code{All.cis.eQTLs}.  The process
of identifying ``best'' cis eQTL per probe leads to a probe-specific
FDR.  In \code{All.cis.eQTLs} we enumerate all probes and all SNP
with FDR at most \code{maxfdr}, not just the best scoring
SNP per probe.}
\item{inbestcis}{Used in \code{All.cis.eQTLs}.  An instance
of \code{\linkS4class{mcwBestCis}} that can be used
to speed up the extraction of All.cis eQTL.}
\item{smFilter4cis}{Used in \code{All.cis.eQTLs}.  
A function accepting and returning an smlSet instance.
When
\code{inbestcis} parameter is NULL, this filter will be used
for identifying the best SNP per probe.}
\item{smFilter4all}{Used in \code{All.cis.eQTLs}.  A function
accepting and returning an smlSet instance.
This filter will be used
for identifying the best SNP per probe.  This filter
should not affect the number of probes.}
  \item{x}{instance of \code{mcwBestCis}}
  \item{type}{character, either 'data.frame' or 'GRanges'}
  \item{excludeRadius}{numeric, defaulting to NULL; if non-null,
defines radius around gene region that is excluded for cis SNP
scoring; must be less than radius}
}
\details{
\code{geneApply} can be set to \code{parallel::mclapply}, for example,
in a multicore context.

\code{mcwBestCis} stands for 'multi-chromosome-wide best cis'
eQTL report container.

It is possible that the filtering processes should
be broken into genotype filtering and expression probe
filtering.

\code{fdr(x)} will return a numeric vector of
plug-in FDR estimates corresponding to probe:association tests
as ordered in the fullreport of a *Cis container.
More metadata should be attached to the output of this
function.
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
