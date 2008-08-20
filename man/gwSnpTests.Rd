\name{gwSnpTests}
\alias{gwSnpTests}
\alias{gwSnpTests,formula,smlSet,cnumOrMissing-method}
\alias{gwSnpTests,formula,smlSet,snpdepth-method}
\title{methods for iterating association tests (expression vs SNP)
   across genomes or chromosomes}
\description{methods for iterating association tests (expression vs SNP)
   across genomes or chromosomes}
\usage{
gwSnpTests(sym, sms, cnum, ...)
}
\arguments{
  \item{sym}{ genesym, probeId, or formula instance}
  \item{sms}{ smlSet instance}
  \item{cnum}{ chrnum instance or missing}
  \item{\dots}{ \dots}
}
\details{
invokes snpMatrix test procedures as appropriate
}
\value{
gwSnpScreenResult or cwSnpScreenResult instance
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
gs1 = GeneSet(c("CPNE1", "RBM12"))
geneIdType(gs1) = SymbolIdentifier()
f2 = gwSnpTests(gs1~male, hmFou, chrnum(20))
f2
names(f2)
plot(f2[["RBM12"]])
# show 'smlSet-wide' fit
f3 = gwSnpTests(gs1~male, hmFou)
f3
}
\keyword{ models }
