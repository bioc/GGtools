\name{GGtools-RangedData}
\docType{methods}
\alias{coerce,cwSnpScreenResult,RangedData-method}
\title{Transform results of gwSnpTests to browser tracks}
\description{Create a browser track from a chromosome-wide SNP screen}
\section{Coercion}{
  \describe{
    \item{}{\code{as(object, "RangedData")}: Coerce a
      \code{\linkS4class{cwSnpScreenResult}}, \code{object}, to a
      \code{\linkS4class{RangedData}} instance, with the genomic
      coordinates -log10 p-values for each SNP
    }
  }
}
\keyword{methods}
