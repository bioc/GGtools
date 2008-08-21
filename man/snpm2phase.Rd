\name{snpm2phase}
\alias{snpm2phase}
\title{convert information in a snp.matrix to PHASE input format  }
\description{convert information in a snp.matrix to PHASE input format
}
\usage{
snpm2phase(snpm, cnum, outfilename)
}
\arguments{
  \item{snpm}{ snp.matrix instance}
  \item{cnum}{ chromosome number as chrnum instance}
  \item{outfilename}{ character name of file to write}
}
\details{ follows phase 2.1 documentation for input format
}
\value{
writes to a file and gives a message
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(smtest)
tt = tempfile()
snpm2phase(smtest, chrnum(20), tt)
readLines(tt)
unlink(tt)
}
\keyword{ models }
