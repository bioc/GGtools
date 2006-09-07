\name{snps}
\alias{snps}
\title{accessor for genotype data in a ggExprSet  }
\description{accessor for genotype data in a ggExprSet
}
\usage{
snps(x)
}
\arguments{
  \item{x}{ instance of \code{\link{racExSet-class}}}
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
data(chr20GGdem)
dim(snps(chr20GGdem))
snps(chr20GGdem)[1:4,1:4]
}
\keyword{ models }
