\name{eqtlTests}
\alias{eqtlTests}
\alias{ieqtlTests}
\alias{meqtlTests}
\alias{cisScores}
\alias{imphm3_1KG_20_mA2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
perform genome x transcriptome eQTL searches with high-performance options
}
\description{
perform genome x transcriptome eQTL searches with high-performance options
}
\usage{
eqtlTests(smlSet, rhs = ~1 - 1, runname = "foo", targdir = "foo", geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, ...)
ieqtlTests(smlSet, rhs = ~1 - 1, rules, runname = "ifoo", targdir = "ifoo", geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, ...)
meqtlTests(listOfSmls, rhslist, runname = "mfoo", targdir = "mfoo", geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, harmonizeSNPs = FALSE, ...)
cisScores( mgr, ffind=1, chr, snpGR, radius=5e5, applier=lapply ) 
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
  \item{\dots}{
parameters passed to \code{\link[snpMatrix]{snp.rhs.tests}}
}
  \item{mgr}{an instance of \code{eqtlTestsManager}}
  \item{ffind}{index into the list of ff files managed by \code{mgr} to be used for
obtaining scores}
  \item{chr}{token identifying chromosome of interest}
  \item{snpGR}{instance of \code{\link[GenomicRanges]{GRanges-class}} with SNP location}
  \item{radius}{numeric: number of basepairs up and downstream of TSS to be used for cis filtering}
  \item{applier}{function, typically either \code{lapply} or \code{mclapply} if multicore
services are desired}
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
imputations to loci not covered in the phase 3 hapmap data in \code{hm3ceuSMS}.
}
\value{
instance of \code{eqtlTestsManager}
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
library(ceuhm3)
if (!exists("hm3ceuSMS")) data(hm3ceuSMS)
library(illuminaHumanv1.db)
cptag = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
indc = which(featureNames(hm3ceuSMS) == cptag[1])
hm = hm3ceuSMS[c(indc,1:19),]  # reduce problem
td = tempdir()
curd = getwd()
setwd(td)
time.lapply = unix.time(e1 <- eqtlTests( hm, ~male ))
e1
dir("foo")
if (require(multicore)) {
 time.mclapply = unix.time(e2 <- eqtlTests(hm, ~male, targdir="bar", runname="bar",
#      geneApply = mclapply ) )
      chromApply = mclapply ) )
 print(time.lapply)
 print(time.mclapply)
}
mm = meqtlTests( list(hm, hm), list(~male, ~male))
mm
mm@fflist[[1]][1:4,1:4]
e1@fflist[[1]][1:4,1:4]  # for nonmissing tests, mm = 2xe1
#e3 = eqtlTests(hm, ~male, targdir="barz", runname="barz",
#      geneApply = mclapply, computeZ = TRUE)
e4 = eqtlTests(hm[1:5,], ~male, runname="abc")
e5 = eqtlTests(hm[6:10,], ~male, runname="def")
#
# create and illustrate use of the director concept 
#  here we take the managers for two gene sets (1:5, 6:10)
#  and create the director that allows queries across them
#
cd = mkCisTransDirector(list(e4,e5), "demo2x5.sqlite", snptabname="snpnames", 
       probetabname="probenames", probeanno=hm@annotation, commonSNPs=TRUE)
#
# first we have specific snp-gene queries
#
cd[c("rs6650104", "rs11510103"), "GI_10047091-S"]
#
# now we want values for two SNP over all genes 
#  there were only a few genes so we can look at all
#
cd[c("rs6650104", "rs11510103"),]
#
# now get all SNP for a given gene -- this is a large request
#
dem2 = cd[, "GI_10047091-S"]
sapply(dem2[[1]], length)
hmi = hm[ chrnum("chr20"), ]
data(imphm3_1KG_20_mA2)
e1i = ieqtlTests(hmi, ~male, rules=imphm3_1KG_20_mA2 )
e1i
topFeats( sym = "CPNE1", mgr=e1i, ffind=1, anno="illuminaHumanv1.db" )
topFeats(probeid="GI_23397697-A", mgr=cd)  
topFeats(rsid="rs6060632", mgr=cd)
#
library(ceu1kg)
data(ceu1kgMeta_20)
library(multicore)
cc = cisScores( e1i, 1, "chr20", ceu1kgMeta_20, applier=mclapply, radius=5e5)
lapply(cc, function(x) sort(x,decreasing=TRUE)[1:3]) # now check scores obtained under permutation
e1ip = ieqtlTests(permEx(hmi), ~male, rules=imphm3_1KG_20_mA2, targdir="p1", runname="p1" )
cc2 = cisScores( e1ip, 1, "chr20", ceu1kgMeta_20, applier=mclapply, radius=5e5)
lapply(cc2, function(x) sort(x,decreasing=TRUE)[1:3])
setwd(curd)
}
\keyword{ models }
