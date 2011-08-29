\name{exome_minp}
\alias{exome_minp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
acquire  minimum p-value for association between genotype
and expression
}
\description{
acquire  minimum p-value for association between genotype
and expression in context of exome genotyping -- where a list
of SNPs associated with genes or exons governs organization of
tests, and minimum p-value per gene or exon is all that is required
}
\usage{
exome_minp(smlSet, fmla, targdir, runname, snpl, feat=NULL, mgr = NULL, scoreApply=lapply, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
basic genotype plus expression structure; this must have an smList() result of length 1
(all SNP in one SnpMatrix regardless of number of chromosomes)
}
  \item{fmla}{
formula expressing covariates to be found in phenoData of smlSet
and used in each association model
}
  \item{targdir}{
folder where ff files will be written
}
  \item{runname}{
prefix for names of ff files
}
  \item{snpl}{
a named list, with one element per gene or exon, each element is name
of snps assayed for the associated gene or exon; names of list elements are the gene
or exon names
}
  \item{feat}{ name of feature for focused reporting; important if names
of features of original smlSet don't agree with names of snpl}
  \item{mgr}{
if an eqtlTestsManager (with fflist of length 1) is already available, this can
be used instead of constructing one from the smlSet
}
  \item{scoreApply}{ lapply-like function to be used to
compute scores -- use mclapply for multicore deployment }
  \item{\dots}{
parameters passed to eqtlTests
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%%  ~Describe the value returned
%%%  If it is a LIST, use
%%%  \item{comp1 }{Description of 'comp1'}
%%%  \item{comp2 }{Description of 'comp2'}
%%% ...
%}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(hmceuB36.2021)
hmlit = hmceuB36.2021[ chrnum(20), ]
library(illuminaHumanv1.db)
cptag = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
indc = which(featureNames(hmlit) == cptag[1])
hm = dropMonomorphies(hmlit[c(indc,1:19),])  # reduce problem
curd = getwd()
td = tempdir()
setwd(td)
sl = colnames(smList(hm)[[1]])[1:80]
sl = split(sl, rep(1:20, each=4))
names(sl) = featureNames(hm)
e1 = exome_minp( hm, ~male, "ex1", "ex1", sl )
e1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
