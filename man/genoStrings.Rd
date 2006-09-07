\name{genoStrings}
\alias{genoStrings}
\title{create a character vector of genotype value strings }
\description{ 
create a character vector of genotype value strings 
}
\usage{
genoStrings(racExSet, rsnum)
}
\arguments{
  \item{racExSet}{ instance of \code{\link{racExSet-class}}}
  \item{rsnum}{ string that is a rowname in snps(racExSet)}
}
\details{
computes unphased genotype code from the rare allele count and other metadata
in the racExSet
}
\value{
character vector
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{For interoperability with the genetics package,
e.g., \code{\link[genetics]{genotype}}  }


%\seealso{  }

\examples{
data(chr20GGdem)
genoStrings(chr20GGdem, "rs6060535")
}
\keyword{ models }
