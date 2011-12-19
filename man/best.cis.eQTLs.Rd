\name{best.cis.eQTLs}
\alias{best.cis.eQTLs}
\alias{fullreport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simple interface to search for cis-eQTL
}
\description{
simple interface to search for cis-eQTL
}
\usage{
best.cis.eQTLs(smpack, fmla, cisRadius = 50000, 
 genome = "hg19", exTransform = function(x) x, 
 folderstem = "cisScratch", nperm = 2, 
 geneApply = lapply, chromApply = lapply, 
 cleanChrn = function(x) gsub("chr", "", x), 
 additionalSNPGR = NULL, useTxDb = FALSE, 
 verbose = TRUE, dropChr = c("X", "Y", "M"))

fullreport(x)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{instance of eqtlFDRSummary, to be described}
  \item{smpack}{
character string naming a package that can be used successfully
with \code{\link[GGBase]{getSS}}: this package has inst/parts
holding \code{\link[snpStats]{SnpMatrix-class}} instances and
data/eset.rda holding an \code{\link[Biobase]{ExpressionSet-class}}
instance named \code{ex}, and is typically created using the
\code{\link{externalize}} function on a unified 
\code{\link[GGBase]{smlSet-class}} instance.
}
  \item{fmla}{
right hand side of a formula; \code{pData} of the \code{smlSet} instances
obtained via \code{getSS} will be used to find variable bindings
}
  \item{cisRadius}{
nucleotide count upstream or downstream of gene coding region within which
SNP associations will be regarded as ``cis''
}
  \item{genome}{
string, either \code{"hg18"} or \code{"hg19"}
}
  \item{exTransform}{
a function that operates on an \code{smlSet} instance and returns one;
will typically operate on the expression data (for example, removing
expression heterogeneity) and/or the genotype data (for example, using
\code{\link{MAFfilter}})
}
  \item{folderstem}{
string where test results will be stored out of memory (in ff archives)
}
  \item{nperm}{
number of permutations of expression against genotype that will be
used for FDR computation
}
  \item{geneApply}{
lapply-like function to be used to iterate over genes within chromosomes
}
  \item{chromApply}{
lapply-like function to be used to iterate over chromosomes
}
  \item{cleanChrn}{
function on string that will harmonize tokens used for referring to
chromosomes -- typically if the \code{smPack} has genotype parts
named \code{chrnn}, we will want to remove the \code{chr} prefix
so that the remainder can be used to look up gene annotation
}
  \item{additionalSNPGR}{
if the SNP annotation resource does not have all SNP locations of interest
(for example, if imputation has been used to generate chrnn:nnnnnn SNP)
this GRanges instance can supply additional SNP identified in \code{smPack}
parts
}
  \item{useTxDb}{
logical variable: if TRUE, gene annotations will be taken from a TxDb.Hsapiens
package corresponding to value of \code{genome}
}
  \item{verbose}{
logical: if TRUE, provide real time details
}
  \item{dropChr}{
string vector identifying parts (such as chromosomes X, Y) that
are not to be analyzed
}
}
\details{
The intention is to simplify cis eQTL searches, relying as much as possible
on Bioconductor annotation resources to determine gene-chromosome mapping and
gene-SNP distances.  Gene and SNP filtering can be accomplished through
the \code{exTransform} parameter.

For each chromosome, the set of probes resident 
thereupon are isolated and all SNP-probe associations are tested using
\code{\link{eqtlTests}}.  This creates a temporary \code{ff} archive with
compactly represented chi-squared statistics for all associations.
These statistics are then harvested using the location information and
radius specification in the call.

Iteration over genes is nested within iteration over chromosomes.  A
potential workflow assigns chromosomes to MPI slaves and then uses
\code{geneApply=mclapply} to allow concurrent computing over genes.

}
\value{
instance of class \code{eqtlFDRSummary}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
This procedure reports on the best cis associated SNP to each
gene.  Other policies for assessing associations can be explored
with \code{policywiseFDRtab}, and, for those genes exhibiting
associations, additional (but not maximally) associated loci may be identified
using \code{gwSnpTests} or \code{eqtlTests}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{eqtlTests}}, \code{\link{gwSnpTests}}, 
\code{\link[GGBase]{smlSet-class}}, \code{\link{externalize}}
}
\examples{
# this is a very small example that will yield 4 hits
\dontrun{
 best.cis.eQTLs("GGdata", ~male, exTransform=function(x)
   nsFilter(clipPCs(MAFfilter(x, lower=.05), 1:10), var.cutoff=.975), 
   dropChr=c("X", "Y", 1:20), folderstem="daga5", geneApply=mclapply)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
