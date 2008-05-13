\name{parse2SM}
\alias{parse2SM}
\title{ parse a set of HapMap downloads into snpMatrix entities }
\description{ parse a set of HapMap downloads into snpMatrix entities }
\usage{
parse2SM(pop="CEU", rel="r23a", bld="b36")
}
\arguments{
  \item{pop}{ string giving the population code used in
the hapmap filename}
  \item{rel}{ string giving minor release info}
  \item{bld}{ string giving build info encoded in hapmap
filenames}
}
\details{
Run in a folder populated with 24 .txt.gz hapmap downloads
for a given population.  Returns a list with a snp.matrix
list and a support (metadata) data.frame, after parsing these
files.
}
\value{
Returns a list with a snp.matrix
list and a support (metadata) data.frame.
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
#
# need a slice to demonstrate
#
}
\keyword{ models }
