\name{makeLocNCDF}
\alias{makeLocNCDF}
\title{support creation of NetCDF representation (external) of snp locations  }
\description{support creation of NetCDF representation (external) of snp locations  }
\usage{
makeLocNCDF(smSupp, fn, report=TRUE, doRSID=FALSE)
}
\arguments{
  \item{smSupp}{  data.frame with SNP support info as returned
by snpMatrix read.HapMap.data facility}
  \item{fn}{ filename to write NetCDF in}
  \item{report}{ should a brief report be written to stdout upon completion?}
  \item{doRSID}{ should dbSNP ID (rs numbers) be included in the store?
defaults to FALSE because they consume a lot of space and are typically
already available in the data under analysis.}
}
\details{
run for the side-effect of creating a NetCDF store.
}
\value{
invisibly returns the result of close.ncdf on the created store.
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
# just infrastructure for developer
}
\keyword{ models }
