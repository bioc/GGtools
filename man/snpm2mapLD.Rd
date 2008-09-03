\name{snpm2mapLD}
\alias{snpm2mapLD}
\title{ prepare input to mapLD function for haplotype block identification }
\description{prepare input to mapLD function for haplotype block identification
}
\usage{
snpm2mapLD(x, chrnum, runMAP=TRUE, ...)
}
\arguments{
  \item{x}{ snp.matrix instance}
  \item{chrnum}{ chromosome number}
  \item{runMAP}{ logical indicating whether or not to run mapLD}
  \item{\dots}{ additional arguments to mapLD}
}
\details{
sets up a data frame suitable for mapLD, and will invoke with
appropriate arguments identifying columns for alleles and
other identifiers if runMAP is TRUE (default).
}
\value{a list with element struc holding the data frame,
and mapLD output if requested.  Note that mapLD writes an eps
file to disk *sigh*.
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(smtest)
ss = snpm2mapLD(smtest, chrnum=20, runMAP=FALSE)
ss
# you could run mapLD on ss[[1]]
}
\keyword{ models }
