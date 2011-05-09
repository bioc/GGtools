\name{cisProxScores}
\alias{cisProxScores}
\alias{scoresByGenes}
%\alias{mcisProxScores}
\alias{interleave2cis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create, combine, and harvest eqtlTestsManager instances
to collect all eQTL tests satisfying certain gene proximity conditions
}
\description{
create, combine, and harvest eqtlTestsManager instances
to collect all eQTL tests satisfying certain gene proximity conditions
}
\usage{
cisProxScores(smlSet, fmla, dradset, direc = NULL, folder, runname, geneApply = lapply, saveDirector = TRUE, snpGRL=NULL,
 geneGRL=NULL, snpannopack="SNPlocs.Hsapiens.dbSNP.20100427", ffind=NULL, ...)

#mcisProxScores (listOfSmlSets, listOfFmlas, dradset, direc = NULL, 
#    folder, runname, geneApply = mclapply, saveDirector = TRUE, 
#    makeCommonSNPs = FALSE, snpGRL=NULL,
#    geneGRL=NULL, snpannopack="SNPlocs.Hsapiens.dbSNP.20100427", ffind=NULL, ...)
#
interleave2cis( cisp, permcisp )

scoresByGenes(cps, intvind = 1, as.GRanges=TRUE, dups2max=TRUE, snpGR=NULL,
  scoreConverter=function(x)x )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
instance of \code{\link[GGBase]{smlSet-class}}
}
  \item{fmla}{
the right-hand side of a standard modeling formula -- no
dependent variable; the expression values in the \code{smlSet}
will be used successively as dependent variables
}
  \item{dradset}{
a numeric vector indicating the boundaries within which
test scores will be tabulated.  For example, if \code{dradset}
is \code{c(5000,10000,25000)} then scores will be tabulated
for SNP in the regions (0-5kb) from start or end of gene,
(5-10kb), (10-25kb).
}
  \item{direc}{
an instance of \code{\link[GGtools]{multiCisDirector-class}}; if non-null,
\code{\link{eqtlTests}} will not be run, but the tests managed
by managers in the \code{direc} instance will be used
}
  \item{folder}{
used to set \code{targdir} parameter when \code{eqtlTests}
is run; ignored if \code{direc} is non-null
}
  \item{runname}{
used to set \code{runname} parameter when \code{eqtlTests}
is run; some mangling will be applied.  Ignored if \code{direc} is non-null
}
  \item{geneApply}{
iteration function (like \code{lapply}) to be used for each
expression probe (gene); passed to \code{eqtlTests}; 
the setting is also used for some annotation-based iterations;
if multicore package 
is present, setting this parameter to \code{mclapply} is advised
}
  \item{saveDirector}{
logical; since it is expensive to compute the \code{multiCisDirector}
that will be harvested, we may want to serialize it; if so set
\code{saveDirector} to TRUE.  If set to true the function stores
an object with name \code{paste(folder,"_director",".rda",sep="")}
in the current working folder.
}
%  \item{geneCentric}{logical to specify whether report is
%oriented towards SNP scores as features of genes, or gene
%scores as features of SNPs.  See details below.}
%  \item{retain}{ numeric used only when \code{geneCentric} is \code{FALSE}.
%tells number of gene scores to retain per SNP.}
  \item{\dots}{
arguments passed to \code{\link{eqtlTests}}
 }
 \item{listOfSmlSets}{for \code{mcisProxScores}, a list of
smlSets that are to be sources for eQTL test scores that will be 
summed}
 \item{listOfFmlas}{for \code{mcisProxScores}, a list of formulas
to be used with \code{snp.rhs.tests}, assumed to be ordered to
correspond to elements of \code{listOfSmlSets}}
 \item{makeCommonSNPs}{for \code{mcisProxScores}, a logical telling
whether the sets of SNPs elements of the \code{listOfSmlSets} should
be reduced to their intersection; this can be slow, and can be
done externally using the function of the same name.}
 \item{snpGRL}{named list of GRanges instances  with SNP locations;
list element names must coincide with names of smList entries in
smlSet}
 \item{geneGRL}{named list of GRanges instances  with gene extents;
list element names must coincide with names of smList entries in
smlSet}
 \item{snpannopack}{string naming package with SNPlocs information}
 \item{cisp}{result of \code{cisProxScores}}
 \item{permcisp}{result of \code{cisProxScores}}
 \item{ffind}{usually 1
  for cis applications where one chromosome of SNP is selected at a time}
 \item{cps}{instance of \code{\link{cisProxScores-class}}
 \item{intvind}{index of cis interval to be evaluated (usually the cisProxScores
has been run with a dradset specifying a set of disjoint intervals, given
by names(cps), for example, where cps is the cisProxScores instance}
 \item{as.GRanges}{logical indicating that scores should be returned bound
to GRanges for SNP addresses}
 \item{dups2max}{logical indicating policy for dealing with SNP that occur
multiple times in this gene-oriented survey -- the duplicate SNP are reduced
to the highest-scoring one}
 \item{snpGR}{GRanges instance with SNP addresses -- need not be identical
to set of SNP analyzed}
 \item{scoreConverter}{a function accepting and returning numeric; for minus log10 p,
use \code{scoreConverter=function(x)-log10(1-pchisq(x,df))}, where df is typically 1}
}  % end arguments{}
\details{
This function computes tests for all same-chromosome eQTL
up to the maximum distance given in \code{dradset}
and returns a named list with chi-squared statistics computed
by \code{\link[snpStats]{snp.rhs.tests}}

The \code{interleave2cis} function helps with general comparison
of distributions of real scores to distributions obtained
after permutation of expression values against genotypes.  See
the example.

%If \code{geneCentric} is TRUE, the logic of reporting is straightforward.
%For each gene, a family of genomic
%intervals relative to the coding interval is specified
%and scores for SNPs are collected for these intervals.
%
%If \code{geneCentric} is FALSE, the family of intervals is used to
%define whether SNP should be included or not in final reporting.
%If a SNP is within the specified radius of some gene, then it is kept.
%All same-chromosome gene scores for each SNP are sorted and the top \code{retain} are
%returned.
}
\value{
a list with one component per `radius' derived from \code{dradset}

each radius-associated component includes a list with
one element per chromosome of the SNP data in the \code{smlSet}

each chromosome-associated sublist includes a list for each
gene mapped to the chromosome, with contents a column-vector
of test results for all SNP within the radius of the enclosing
component; see the example for further concreteness 
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{eqtlTests}}
}
\examples{
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
hm = hmceuB36.2021
td = tempdir()
cd = getwd()
on.exit(setwd(cd))
setwd(td)
library(illuminaHumanv1.db)
g20 = intersect(get("20", revmap(illuminaHumanv1CHR)),
   featureNames(hm))[1:10]
g21 = intersect(get("21", revmap(illuminaHumanv1CHR)),
   featureNames(hm))[1:10]
hm = hm[probeId(c(g20,g21)), ] # restrict to small number of genes
try(unlink("man", recursive=TRUE))  # in tempdir
set.seed(1234)  # necessary for dealing with null imputation of missing
f1 = cisProxScores( hm, ~male, c(5000,10000,25000), folder="man",
   runname="man", geneApply=lapply, ffind=1 )
length(f1)  # number of proximity regions specified in dradset
length(f1[[1]]) # number of chromosomes of SNP data in smlSet
length(f1[[1]][[1]]) # number of genes in smlSet
                     # mapping to first chromosome in smlSet 
                     # SNP data
length(f1[[1]][[2]]) # number of genes mapping to second chr...
sapply(f1, function(x)max(unlist(x)))
sapply(f1, function(x)length(x[[1]]))
lapply(f1, function(x)names(x[[1]]))
lapply(f1, function(x)rownames(x[[1]][[1]][[1]]))
set.seed(1234)
try(unlink("pman", recursive=TRUE))  # in tempdir
pf1 = cisProxScores( permEx(hm), ~male, c(5000, 10000, 25000), folder="pman",
  runname="pman", geneApply=lapply, ffind=1)
i1o = interleave2cis( f1, pf1 )
opar = par(no.readonly=TRUE)
par(las=2, mar=c(12, 5, 5, 5))
boxplot(lapply(i1o, unlist), range=0, main="compare observed to expr-permuted eQTL test scores")
par(opar)
load("man_director.rda")
man_director
\dontrun{
set.seed(1234)  # necessary for dealing with null imputation of missing
mm = mcisProxScores( list(hm,hm), list(~male,~male), 
  dradset=c(5000,10000,25000), folder="mmm", runname="MMM", ffind=1)
}

setwd(cd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
