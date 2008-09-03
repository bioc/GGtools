\name{invokePhase-methods}
\docType{methods}
\alias{invokePhase}
\alias{invokePhase-methods}
\alias{invokePhase,phaseInput,chrnum,character,character,character,logical-method}
\alias{invokePhase,snp.matrix,chrnum,character,character,character,logical-method}
\title{ ~~ Methods for Function invokePhase in Package `GGtools' ~~}
\description{
 ~~ Methods for function \code{invokePhase} in Package `GGtools' ~~
}
\section{Methods}{
\describe{
\item{x = "snp.matrix", cnum = "chrnum", parmstring = "character", globpname = "character", where2run = "character", doParse = "logical"}{ 
  transform snp.matrix entity to phaseInput (uses tempfile()) and
   invokes PHASE
  }

\item{x = "phaseInput", cnum = "chrnum", parmstring = "character", globpname = "character", where2run = "character", doParse = "logical"}{ 
  for prepared 'phaseInput' structure, invoke PHASE
  }

}}
\examples{
data(smtest)
invokePhase(smtest, chrnum(20), "", 
  "/Users/stvjc/Desktop/CompCVS/PHASE/MAC/PHASE", ".", TRUE)
}
\keyword{methods}
