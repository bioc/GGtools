\name{transScores}
\alias{transScores}
\alias{transScores.legacy}
\alias{meta.transScores}
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
transScores( tconfig )

transScores.legacy(smpack, snpchr = "chr1", rhs, K = 20, targdirpref = "tsco", geneApply = lapply, 
  chrnames = paste("chr", as.character(1:22), sep = ""), geneRanges = NULL, snpRanges = NULL, 
  radius = 2e+06, renameChrs = NULL, probesToKeep = NULL, batchsize = 200, 
  genegran = 50, shortfac = 10, wrapperEndo = NULL, 
  geneannopk = "illuminaHumanv1.db", 
  snpannopk = snplocsDefault(), gchrpref = "",
                 schrpref = "ch", exFilter=function(x)x,
   smFilter=function(x)x,
   SSgen=GGBase::getSS)

meta.transScores (smpackvec = c("GGdata", "hmyriB36"),
    snpchr = "22", rhsList=list(~1, ~1), K = 20, targdirpref = "mtsco",
    geneApply = lapply, chrnames = as.character(21:22),
    radius = 2e+06, renameChrs=NULL,
    probesToKeep=NULL, batchsize=200, genegran=50, shortfac=10, wrapperEndo=NULL,
    geneannopk = "illuminaHumanv1.db", snpannopk = snplocsDefault(),
    gchrpref = "", schrpref="ch",
    exFilterList= list(function(x)x, function(x)x),
    SMFilterList = list(function(x)x, function(x)x),
    SSgen = GGBase::getSS)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tconfig}{instance of \code{\link{TransConfig-class}}}
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
  \item{rhsList}{
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
passed to \code{\link[GGBase:make_smlSet]{getSS}}
}
  \item{probesToKeep}{
passed to \code{\link[GGBase:make_smlSet]{getSS}}
}
  \item{batchsize}{
defines batch size for \code{\link[ff:ffapply]{ffrowapply}}
}
  \item{genegran}{
passed to \code{\link{eqtlTests}}
}
  \item{shortfac}{
passed to \code{\link{eqtlTests}}
}
  \item{wrapperEndo}{
a function accepting and returning an smlSet instance, evaluated
before numerical analysis of associations, which will be executed
on the output of this function
}
  \item{gchrpref}{
prefix to convert \code{chrnames} into appropriate tokens for
obtaining gene metadata; in future this may need to be a string
transformation function
}
  \item{schrpref}{
prefix to convert \code{chrnames} into appropriate tokens for
use with \code{snplocs} for the SNP location information
package identified in \code{snpannopack} parameter below
}
  \item{geneannopk}{
character string naming a Bioconductor .db expression chip annotation package
}
  \item{snpannopk}{
character string naming a Bioconductor SNPlocs.* SNP metadata package
}
  \item{exFilter}{function to transform ExpressionSet component of
retrieved smlSet, to reduce probe sets in use, for example}
  \item{smFilter}{function to transform smlSet instance before use; filter can
affect genotypes in smList(x)[[1]], for example}
    \item{exFilterList}{list of functions serving as exFilters for
each of the elements of smpackvec}
    \item{SMFilterList}{list of functions servicing as wrapperEndos
for each of the elements of smpackvec}
 \item{SSgen}{function to be used to create smlSet instance
for testing -- in general, GGBase::getSS has been used to pull
the ExpressionSet and SnpMatrix data from a named package, but
in some cases a specialize task is needed to create the desired
smlSet.  Whatever is passed to SSgen must return an smlSet instance.}
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
