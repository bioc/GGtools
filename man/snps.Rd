\name{snps}
\alias{snps}
\title{accessor for genotype data in a ggExprSet  }
\description{accessor for genotype data in a ggExprSet
}
\usage{
snps(x)
}
\arguments{
  \item{x}{ instance of \code{\link{ggExprSet}}}
}
\details{
simple accessor, can be inefficicent
}
\value{
instance of data.frame
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(c20GGceu)
dim(snps(c20GGceu))
snps(c20GGceu)[1:4,1:4]
}
\keyword{ models }
