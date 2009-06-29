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
  \item{sms}{ \link[GGBase]{smlSet} instance}
  \item{cnum}{ chrnum instance or missing}
  \item{cs}{ chunksize specification }
  \item{\dots}{ \dots}
}
\details{
invokes \code{snpMatrix} package test procedures 
(e.g., \code{\link[snpMatrix]{snp.rhs.tests}} as appropriate


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
\code{\link[GGBase]{cwSnpScreenResult-class}} instance
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
plot(f1)
# show how to avoid adjusted fit
f1b = gwSnpTests(genesym("CPNE1")~1-1, hmFou, chrnum(20))
# show gene set modeling on chromosome
library(GSEABase)
gs1 = GeneSet(c("CPNE1", "ADA"))
geneIdType(gs1) = SymbolIdentifier()
f2 = gwSnpTests(gs1~male, hmFou, chrnum(20))
f2
names(f2)
plot(f2[["ADA"]])
# show 'smlSet-wide' fit
f3 = gwSnpTests(gs1~male, hmFou)
f3
# now use a phenoVar
f3b = gwSnpTests(phenoVar("persid")~male, hmFou, chrnum(20))
topSnps(f3b)
\dontrun{
# in example() we run into a problem with sys.call(2); works
# in interpreter
f4 = gwSnpTests(gs1~male, hmFou, snpdepth(250), chunksize(1))
f4
  }
}
\keyword{ models }
