\name{gwSnpTests}
\alias{gwSnpTests}
\alias{residTests}
\alias{residTests,cwSnpScreenResult,smlSet,formula,missing-method}
\alias{chunksize}
\alias{chunksize-class}
\alias{gwSnpTests,formula,smlSet,cnumOrMissing-method}
\alias{gwSnpTests,formula,smlSet,cnumOrMissing,missing-method}
\alias{gwSnpTests,formula,smlSet,cnumOrMissing,ANY-method}
\alias{gwSnpTests,formula,smlSet,snpdepth-method}
\alias{gwSnpTests,formula,smlSet,snpdepth,missing-method}
\alias{gwSnpTests,formula,smlSet,snpdepth,ANY-method}
\alias{gwSnpTests,formula,smlSet,snpdepth,chunksize-method}
\title{methods for iterating association tests (expression vs SNP)
   across genomes or chromosomes}
\description{methods for iterating association tests (expression vs SNP)
   across genomes or chromosomes}
\usage{
gwSnpTests(sym, sms, cnum, cs, ...)
}
\arguments{
  \item{sym}{ genesym, probeId, or formula instance}
  \item{sms}{ \link[GGBase:smlSet-class]{smlSet} instance}
  \item{cnum}{ chrnum instance or missing}
  \item{cs}{ chunksize specification }
  \item{\dots}{ \dots}
}
\details{
invokes \code{snpStats} package test procedures 
(e.g., \code{\link[snpStats]{snp.rhs.tests}} as appropriate


\code{chunksize} can be specified to divide
task up into chunks of chromosomes; \code{gc()} will be run
between each chunk -- this may lead to some benefits
when memory capacity is exceeded

The dependent variable in the formula can have class
genesym (chip annotation package used for lookup),
probeId (direct specification using chip annotation vocabulary),
or phenoVar (here we use a phenoData variable as dependent variable).
If you want to put expression values on the right-hand side
of the model, add them to the phenoData and enter them in the
formula.
}
\value{
\code{\link[GGBase]{gwSnpScreenResult-class}} or 
\code{\link[GGBase:gwSnpScreenResult-class]{cwSnpScreenResult-class}} instance
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
# condense to founders only
hmFou = hmceuB36.2021[, which(hmceuB36.2021$isFounder)]
# show basic formula fit
f1 = gwSnpTests(genesym("CPNE1")~male, hmFou, chrnum(20))
f1
#The following code will create a view of the UCSC
#genome browser:
#if (interactive()) {
#library(rtracklayer)
#f1d = as(f1, "RangedData")
#s1 = browserSession("UCSC")
#s1[["CPNE1"]] = f1d
#v1 = browserView(s1, GenomicRanges(30e6, 40e6, "chr20"), full="CPNE1")
#}
# R-based visualization
#plot(f1) -- no longer supported, need to supply location data -- consider eqtlTests/manhPlot
# show how to avoid adjusted fit
f1b = gwSnpTests(genesym("CPNE1")~1-1, hmFou, chrnum(20))
# show gene set modeling on chromosome
\dontrun{
library(GSEABase)  # functionality abandoned
gs1 = GeneSet(c("CPNE1", "ADA"))
geneIdType(gs1) = SymbolIdentifier()
f2 = gwSnpTests(gs1~male, hmFou, chrnum(20))
f2
names(f2)
#plot(f2[["ADA"]])
# show 'smlSet-wide' fit
f3 = gwSnpTests(gs1~male, hmFou)
f3
}
# now use a phenoVar
f3b = gwSnpTests(phenoVar("persid")~male, hmFou, chrnum(20))
topSnps(f3b)
\dontrun{
# in example() we run into a problem with sys.call(2); works
# in interpreter
f4 = gwSnpTests(gs1~male, hmFou, snpdepth(250), chunksize(1))
f4
#
  }
# illustrate alternate approach to expression feature enumeration
#
\dontrun{  # nice but out of scope
data(smlSet.example)
esml = as(smlSet.example, "ExpressionSet")
library(genefilter)
annotation(esml) = "illuminaHumanv1" # drop .db
library(illuminaHumanv1.db)
fesml = nsFilter(esml)[[1]] # unique entrez ids + other filters
fn = featureNames(fesml)
eids = unlist(mget(fn, illuminaHumanv1ENTREZID))
featureNames(fesml) = as.character(eids)
fesml = make_smlSet( fesml, smList(smlSet.example) )
# now we have an smlSet with Entrez ID featureNames
annotation(fesml) = "org.Hs.eg"
mygs = GeneSet(c("ZNF253", "MRS2"), geneIdType = SymbolIdentifier())
geneIdType(mygs) = AnnotationIdentifier("org.Hs.eg")
tt = gwSnpTests(mygs~male, fesml)
lapply(tt, topSnps)
}
}
\keyword{ models }
