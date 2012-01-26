\name{eqtlTests}
\alias{eqtlTests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute association statistics between all probes and SNP in an
smlSet instance
}
\description{
compute association statistics between all probes and SNP in an
smlSet instance, using out-of-memory storage
}
\usage{
eqtlTests(smlSet, rhs = ~1 - 1, runname = "foo", 
	targdir = "foo", geneApply = lapply, 
	shortfac = 100, 
	checkValid = TRUE, useUncertain = TRUE, 
	glmfamily = "gaussian")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
instance of \code{\linkS4class{smlSet}}
}
  \item{rhs}{
fragment of a standard formula, minus a dependent variable (i.e., starts
with tilde); bindings will be sought in \code{pData(smlSet)}
}
  \item{runname}{
string used to identify output ff files
}
  \item{targdir}{
string naming the folder where ff outputs will reside
}
  \item{geneApply}{
analog to \code{lapply} to drive iteration over probes
}
  \item{shortfac}{
ff contents will be multiplied by this quantity and stored as
short integers
}
  \item{checkValid}{
logical, will apply validObject to \code{smlSet} if TRUE
}
  \item{useUncertain}{
logical, passed as \code{uncertain} parameter to 
\code{\link[snpStats]{snp.rhs.tests}} to specify whether uncertain
genotypes will be used (as 'dosage' in GLM fitting)
}
  \item{glmfamily}{
family specification for \code{\link[snpStats]{snp.rhs.tests}}
}
}
\details{
The purpose of this function is to allow very substantial eQTL search 
processes to occur with R.  For several million SNP and tens of thousands
of probes, the storage of test results requires attention to parsimony.
The storage occurs out of memory, using the ff package, and employs
short integers to represent chi squared statistics.  These are scaled
up prior to storage, and will be scaled down prior to use.
}
\value{
returns an instance of \code{eqtlTestsManager}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
hm2ceuSMS = getSS("GGtools", c("20"), renameChrs=c("chr20"))
library(illuminaHumanv1.db)
cptag = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
indc = which(featureNames(hm2ceuSMS) == cptag[1])
#
# get a set of additional genes on chr20
all20 = get("20", revmap(illuminaHumanv1CHR))
g20 = unique(c(all20[1:10], cptag))
#
hm = hm2ceuSMS[probeId(g20),]  # reduce problem
td = tempdir()
curd = getwd()
setwd(td)
time.lapply = unix.time(e1 <- eqtlTests( hm, ~male ))
time.lapply
e1
# best chisq(1) for CPNE1
topFeats(probeId(cptag), e1)
setwd(curd)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
