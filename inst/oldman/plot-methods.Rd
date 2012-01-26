\name{plot-methods}
\docType{methods}
\alias{plot,cwSnpScreenResult,missing-method}
\alias{plot,filteredGwSnpScreenResult,ANY-method}
\alias{plot,filteredMultiGwSnpScreenResult,ANY-method}
\alias{plot,snp.reg.imputation,missing-method}
\title{  Methods for Function plot in Package `GGtools' }
\description{
  Methods for function \code{plot} in Package `GGtools' 
}
\section{Methods}{
\describe{

\item{x = "cwSnpScreenResult", y = "missing"}{ shows results of chromosome-wide screen for expression-associated SNP assuming locations from SNPlocs.Hsapiens.dbSNP.20090506}

\item{x = "cwSnpScreenResult", y = "character"}{ shows results of chromosome-wide screen for expression-associated SNP assuming locations from SNPlocs package named in 'y'}

\item{x = "filteredGwSnpScreenResult", y = "ANY"}{ shows results of genome-wide screen for expression-associated SNP}

\item{x = "filteredMultiGwSnpScreenResult", y = "ANY"}{ fails, need to pick gene at this time }

}}
\keyword{methods}
