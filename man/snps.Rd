\name{snps}
\alias{snps}
\alias{snps,racExSet,missing-method}
\alias{snps,smlSet,chrnum-method}
\title{accessor for genotype data in a ggExprSet  }
\description{accessor for genotype data in a ggExprSet
}
\usage{
snps(x,chr)
}
\arguments{
  \item{x}{ instance of \code{\link{racExSet-class}} or
of \code{\link{smlSet-class}}}
  \item{chr}{ instance of \code{\link{chrnum-class}}}
}
\details{
simple accessor, can be inefficient
}
\value{
instance of data.frame
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(chr20GGdem)
dim(snps(chr20GGdem))
snps(chr20GGdem)[1:4,1:4]
}
\keyword{ models }
