\name{eqtlTests}
\alias{eqtlTests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
perform genome x transcriptome eQTL searches with high-performance options
}
\description{
perform genome x transcriptome eQTL searches with high-performance options
}
\usage{
eqtlTests(smlSet, rhs = ~1 - 1, runname = "foo", targdir = "foo", geneApply = lapply, chromApply = lapply, shortfac = 100, computeZ = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
instance of \code{\link[GGBase]{smlSet-class}}
}
  \item{rhs}{
standard formula without dependent variable; predictors must be found in \code{pData(smlSet)}
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
  \item{\dots}{
parameters passed to \code{\link[snpMatrix]{snp.rhs.tests}}
}
}
\details{
\code{\link[snpMatrix]{snp.rhs.tests}} is run for all genes enumerated in
\code{featureNames(smlSet)} individually as dependent variables, and all SNP  in \code{smList(smlSet)}
as predictors, one by one.  Each model fitted for SNP genotype is additionally
adjusted for elements in \code{rhs}.  There are consequently \code{G*S} test results
where \code{G} is the number of features in \code{exprs(smlSet)}, and \code{S} is the
total number of SNP in \code{smlSet}.  These are stored in \code{ff} files in folder
\code{targdir}.
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
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
library(illuminaHumanv1.db)
cptag = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
indc = which(featureNames(hmceuB36.2021) == cptag[1])
hm = hmceuB36.2021[c(indc,1:19),]  # reduce problem
td = tempdir()
curd = getwd()
setwd(td)
time.lapply = unix.time(e1 <- eqtlTests( hm, ~male ))
e1
dir("foo")
if (require(multicore)) {
 time.mclapply = unix.time(e2 <- eqtlTests(hm, ~male, targdir="bar", runname="bar",
      geneApply = mclapply ) )
 print(time.lapply)
 print(time.mclapply)
}
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
cd[c("rs885550", "rs169757"), "GI_10047099-S"]
#
# now we want values for two SNP over all genes 
#  there were only a few genes so we can look at all
#
cd[c("rs885550", "rs169757"),]
#
# now get all SNP for a given gene -- this is a large request
#
dem2 = cd[, "GI_10047099-S"]
sapply(dem2[[1]], length)
setwd(curd)
}
\keyword{ models }
