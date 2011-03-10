\name{eqtlTests}
\alias{eqtlTests}
\alias{ieqtlTests}
\alias{meqtlTests}
\alias{meqtlTests2}
%\alias{cisScores}
%\alias{imphm3_1KG_20_mA2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
perform genome x transcriptome eQTL searches with high-performance options
}
\description{
perform genome x transcriptome eQTL searches with high-performance options
}
\usage{
eqtlTests(smlSet, rhs = ~1 - 1, runname = "foo", targdir = "foo", 
	geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, 
        checkValid = TRUE, saveSummaries = TRUE, uncert=TRUE, family, genegran=50, ...)
ieqtlTests(smlSet, rhs = ~1 - 1, rules, runname = "ifoo", targdir = "ifoo", 
	geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, uncert=TRUE, saveSummaries=TRUE,
        family, ...)
meqtlTests(listOfSmls, rhslist, runname = "mfoo", targdir = "mfoo", 
	geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, harmonizeSNPs = FALSE, uncert=TRUE, saveSummaries=TRUE, family, ...)
#cisScores( mgr, ffind=1, chr, snpGR, radius=5e5, applier=lapply, minMAF = 0, minGTF = 0, useEnd=FALSE ) 
meqtlTests2(obpaths, rhslist, runname = "mfoo", targdir = "mfoo", 
	geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, harmonizeSNPs = FALSE, uncert=TRUE, saveSummaries=TRUE, family, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
instance of \code{\link[GGBase]{smlSet-class}}
}
  \item{listOfSmls}{
list of instances of \code{\link[GGBase]{smlSet-class}}
}
  \item{rhs}{
standard formula without dependent variable; predictors must be found in \code{pData(smlSet)}
}
  \item{rhslist}{
List of standard formulae without dependent variable; predictors for each formula must be found in \code{pData(smlSet)} associated
with each rhslist element.  In other words, the pData of the kth smlSet in listOfSmls provides the
data to resolve the symbols in the kth formula of rhslist
}
  \item{rules}{instance of \code{\link[snpStats]{ImputationRules-class}} expressing
rules by which unobserved SNP are `imputed' (that is, the value used is the
conditional expectation of
B copy number, which is real-valued and may lie outside [0,2])
}
  \item{runname}{
arbitrary character string that will identify a serialized object storing references to
results
}
  \item{targdir}{
arbitrary character string that will name a folder where results are stored as \code{\link[ff]{ff}}
files
}
  \item{geneApply}{
\code{lapply}-like function for iterating over genes
}
  \item{chromApply}{
\code{lapply}-like function for iterating over chromosomes
}
  \item{shortfac}{
quantity by which chisquared tests will be inflated before coercion to short int
}
  \item{computeZ}{
logical to direct calculation of Zscore instead of X2
}
  \item{harmonizeSNPs}{
logical: it can be time consuming to harmonize SNPs across a long listOfSmls, so you can
do this outside of the function and set harmonizeSNPs=FALSE; if TRUE, it will be done
before statistical processing of the data in this function.
}
  \item{checkValid}{logical: shall the function run validObject on input smlSet?}
  \item{saveSummaries}{logical: shall a set of ff files be stored that includes genotype
and allele frequency data for downstream filtering?}
  \item{uncert}{setting for value of \code{uncertain} argument in \code{\link[snpStats]{snp.rhs.tests}}}
  \item{family}{specify the GLM family to use; defaults to 'gaussian' if left missing}
  \item{\dots}{
parameters passed to \code{\link[snpStats]{snp.rhs.tests}}}
 \item{obpaths}{vector of strings giving paths to serialized smlSet instances to avoid excess consumption of
RAM.  These will be loaded as needed.  }
  \item{mgr}{an instance of \code{eqtlTestsManager}}
  \item{ffind}{index into the list of ff files managed by \code{mgr} to be used for
obtaining scores}
  \item{chr}{token identifying chromosome of interest}
  \item{snpGR}{instance of \code{\link[GenomicRanges]{GRanges-class}} with SNP location}
  \item{radius}{numeric: number of basepairs up and downstream of TSS to be used for cis filtering}
  \item{minMAF}{numeric: minimum minor allele frequency for inclusion of SNP in reporting}
  \item{minGTF}{numeric: minimum genotype frequency for inclusion of SNP in reporting}
  \item{applier}{function, typically either \code{lapply} or \code{mclapply} if multicore
services are desired}
\item{useEnd}{logical to request cis relative to radius
past entire gene sequence
as opposed to local to TSS}
 \item{genegran}{numeric value of frequency at which gene names will be catted to stdout in case
  \code{options()$verbose == TRUE}}
}
\details{
\code{\link[snpStats]{snp.rhs.tests}} is run for all genes enumerated in
\code{featureNames(smlSet)} individually as dependent variables, and all SNP  in \code{smList(smlSet)}
as predictors, one by one.  Each model fitted for SNP genotype is additionally
adjusted for elements in \code{rhs}.  There are consequently \code{G*S} test results
where \code{G} is the number of features in \code{exprs(smlSet)}, and \code{S} is the
total number of SNP in \code{smlSet}.  These are stored in \code{ff} files in folder
\code{targdir}.

\code{imphm3_1KG_20_mA2} is a set of imputation rules for SNP on chromosome 20, where the
1000 genomes genotypes distributed in `pilot1' VCF files are used to create
imputations to loci not covered in the phase 3 hapmap data in \code{hm3ceuSMS}.

cisScores will fail if genes are present that are not on the chromosome
for which scores are requested.
}
\value{
(i,m)eqtlTests returns instance of \code{eqtlTestsManager}

cisScores returns list with elements for each gene consisting of
chi-squared statistics for SNP cis to the genes according to settings
of radius and useEnd
}
%\references{
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
We are using \code{\link{ff}} to manage the extremely voluminous results
of comprehensive eqtl searches with one short int per test.  We do not have 
an approach to handling NA in this framework, so for any nonexistent test
result (due for example to monomorphy or total missingness) we impute a value
from the null distribution of the test statistic being computed -- chisq of
one d.f..  There is no practical risk of misinterpreting such results in
contexts of interest, but this saves us the complication of dealing with
artificial masses of test statistic distributions at zero, for example.

The \code{\link{topFeats}} methods have \code{minMAF} and \code{minGTF}
parameters to assist in filtering results to SNPs
with certain properties; the metadata used for these is
stored in a summary ff structure.  When \code{meqtlTests} is used,
these summary statistics are generated based on only the first
smlSet in the list passed into meqtlTests.  If significant
use cases for these filtering parameters emerge, the summary
statistics will have to be specialized, but listwise summarization is
not really relevant as with meqtlTests we are marginalizing over
populations.  Some projection of the list-wise summaries will be most 
convenient to implement.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(ceuhm3)
if (!exists("hm3ceuSMS")) data(hm3ceuSMS)
library(illuminaHumanv1.db)
cptag = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
indc = which(featureNames(hm3ceuSMS) == cptag[1])
#
# get a set of additional genes on chr20
all20 = get("20", revmap(illuminaHumanv1CHR))
g20 = unique(c(all20[1:10], cptag))
#
hm = hm3ceuSMS[probeId(g20),]  # reduce problem
hm = hm[ chrnum(c("chr20", "chr21")), ]
td = tempdir()
curd = getwd()
setwd(td)
time.lapply = unix.time(e1 <- eqtlTests( hm, ~male ))
e1
dir("foo")
setwd(curd)
#
# additional examples are in the 'extras' folder, extrExt.txt
# 
}
\keyword{ models }
