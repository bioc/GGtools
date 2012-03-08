\name{transScores}
\alias{transScores}
\alias{mtransScores}
\alias{tr1_obs}
\alias{tr1_perm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
obtain the top trans associations for each SNP in an smlSet
}
\description{
obtain the top trans associations for each SNP in an smlSet
}
\usage{
transScores(smpack, snpchr = "chr1", rhs, K = 20, targdirpref = "tsco", geneApply = lapply, 
  chrnames = paste("chr", as.character(1:22), sep = ""), geneRanges = NULL, snpRanges = NULL, 
  radius = 2e+06, renameChrs = NULL, probesToKeep = NULL, batchsize = 200, 
  genegran = 50, shortfac = 10, wrapperEndo = NULL, 
  geneannopk = "illuminaHumanv1.db", 
  snpannopk = "SNPlocs.Hsapiens.dbSNP.20110815", gchrpref = "",
                 schrpref = "ch")

mtransScores (smpackvec, snpchr = "chr1", rhslist, K = 20, targdirpref = "multtsco",
    geneApply = lapply, chrnames = paste("chr", as.character(1:22), sep=""),
    geneRanges = NULL, snpRanges = NULL, radius = 2e+06, renameChrs=NULL,
    batchsize=200, genegran=50, probesToKeep=NULL, shortfac=10, wrapperEndo=NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smpack}{
name of package holding eset.rda providing 'ex' ExpressionSet when loaded, and holding
SnpMatrix instances in inst/parts
}
  \item{smpackvec}{
vector of names of package holding eset.rda providing 'ex' ExpressionSet when loaded, and holding
SnpMatrix instances in inst/parts
}
  \item{snpchr}{
name or vector of chromosome names of SNPs of interest
}
  \item{rhs}{
right hand side of snp.rhs.tests model for which expression is left hand side, e.g., covariates other than genotype
}
  \item{rhslist}{
list of right hand side of snp.rhs.tests model for which expression is left hand side, e.g., covariates other than genotype,
one per element of \code{smpackvec}
}
  \item{K}{
number of most highly associated features to be retained
}
  \item{targdirpref}{
prefix of target folder name (passed to \code{\link{eqtlTests}}
}
  \item{geneApply}{
passed to \code{\link{eqtlTests}}
}
  \item{chrnames}{
names of chromosomes harboring genes that will be tested for association with genotype
}
  \item{geneRanges}{
list of \code{\link[GenomicRanges]{GRanges-class}} instances containing chromosomal coordinate defined regions occupied by genes, with regions partitioned by
chromosomes, and list element names as given in \code{chrnames} above
}
  \item{snpRanges}{
list of \code{\link[GenomicRanges]{GRanges-class}} instances with SNP addresses
}
  \item{radius}{
radius within which an association is considered cis and therefore the corresponding
test statistic is set to zero
}
  \item{renameChrs}{
passed to \code{\link[GGBase]{getSS}}
}
  \item{probesToKeep}{
passed to \code{\link[GGBase]{getSS}}
}
  \item{batchsize}{
defines batch size for \code{\link[ff]{ffrowapply}}
}
  \item{genegran}{
passed to \code{\link{eqtlTests}}
}
  \item{shortfac}{
passed to \code{\link{eqtlTests}}
}
  \item{wrapperEndo}{
a function accepting and returning an smlSet instance
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
  \item{geneannopk}{
character string naming a Bioconductor .db expression chip annotation package
}
  \item{snpannopk}{
character string naming a Bioconductor SNPlocs.* SNP metadata package
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a list with elements
\item{scores}{an S by K ff matrix where S is number of SNPs, K is number of best features to
be retained, with element s,k the kth largest score statistic among association tests
computed for SNP s}
\item{inds}{an S by K ff matrix with s,k element telling which element of guniv (see below) is the gene
giving the kth largest score statistic for association}
\item{guniv}{the vector of gene identifiers defining the universe of genes tested}
\item{snpnames}{vector of SNP identifiers}
\item{call}{the call used to create the result}

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
library(GGdata)
# need to define the geneRanges and snpRanges ...
transScores("GGdata", "20", renameChrs="chr20", chrnames="chr21")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
