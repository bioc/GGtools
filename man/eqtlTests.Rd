\name{eqtlTests}
\alias{eqtlTests}
\alias{ieqtlTests}
\alias{meqtlTests}
\alias{cisScores}
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
        checkValid = TRUE, saveSummaries = TRUE, uncert=TRUE, family, ...)
ieqtlTests(smlSet, rhs = ~1 - 1, rules, runname = "ifoo", targdir = "ifoo", 
	geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, uncert=TRUE, family, ...)
meqtlTests(listOfSmls, rhslist, runname = "mfoo", targdir = "mfoo", 
	geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, harmonizeSNPs = FALSE, uncert=TRUE, family, ...)
cisScores( mgr, ffind=1, chr, snpGR, radius=5e5, applier=lapply, minMAF = 0, minGTF = 0, useEnd=FALSE ) 
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
  \item{rules}{instance of \code{\link[snpMatrix]{snp.reg.imputation-class}} expressing
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
  \item{uncert}{setting for value of \code{uncertain} argument in \code{\link[snpMatrix]{snp.rhs.tests}}}
  \item{family}{specify the GLM family to use; defaults to 'gaussian' if left missing}
  \item{\dots}{
parameters passed to \code{\link[snpMatrix]{snp.rhs.tests}}
}
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
}
\details{
\code{\link[snpMatrix]{snp.rhs.tests}} is run for all genes enumerated in
\code{featureNames(smlSet)} individually as dependent variables, and all SNP  in \code{smList(smlSet)}
as predictors, one by one.  Each model fitted for SNP genotype is additionally
adjusted for elements in \code{rhs}.  There are consequently \code{G*S} test results
where \code{G} is the number of features in \code{exprs(smlSet)}, and \code{S} is the
total number of SNP in \code{smlSet}.  These are stored in \code{ff} files in folder
\code{targdir}.

\code{imphm3_1KG_20_mA2} is a set of imputation rules for SNP on chromosome 20, where the
1000 genomes genotypes distributed in `pilot1' VCF files are used to create
imputations to loci not covered in the phase 3 hapmap data in \code{hm3ceuSMS}.  Not useful
in Bioconductor 2.7 as snpMatrix package lost its imputation facilities just prior to release;
use snpMatrix2 in Bioconductor 2.8.

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
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
library(illuminaHumanv1.db)
cptag = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
indc = which(featureNames(hmceuB36.2021) == cptag[1])
#
# get a set of additional genes on chr20
all20 = get("20", revmap(illuminaHumanv1CHR))
g20 = unique(c(all20[1:10], cptag))
#
hm = hmceuB36.2021[probeId(g20),]  # reduce problem
hm = hm[chrnum(c("20", "21")),]
td = tempdir()
curd = getwd()
setwd(td)
time.lapply = unix.time(e1 <- eqtlTests( hm, ~male ))
e1
dir("foo")
#
# now show how to extract scores cis to genes
#
 if (require(SNPlocs.Hsapiens.dbSNP.20100427)) {
 c20 = getSNPlocs("ch20")
 library(GenomicRanges)
 c20r = GRanges(seqnames="chr20", IRanges(c20$loc, width=1))
 names(c20r) = paste("rs", c20$RefSNP_id, sep="")
 chkcc = cisScores(e1, 1, "chr20", c20r, radius=2e6) 
 sapply(chkcc, max)
}
setwd(curd)
#
# additional examples are in the 'extras' folder, extrExt.txt
# 
}
\keyword{ models }
