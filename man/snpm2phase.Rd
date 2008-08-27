\name{snpm2phase}
\alias{snpm2phase}
\alias{parsePhPairs}
\alias{personalHap}
\title{convert information in a snp.matrix to PHASE input format  }
\description{convert information in a snp.matrix to PHASE input format
}
\usage{
snpm2phase(snpm, cnum, outfilename)
parsePhPairs(fn)
personalHap(x)
}
\arguments{
  \item{snpm}{ snp.matrix instance}
  \item{cnum}{ chromosome number as chrnum instance}
  \item{outfilename}{ character name of file to write}
  \item{fn}{ character name of PHASE pairs file to read}
  \item{x}{ output of parsePhPairs}
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
pin = snpm2phase(smtest, chrnum(20), tt)
pin
readLines(tt)
unlink(tt)
pp = parsePhPairs(system.file("phaseOut/cpne1_20k.out_pairs", package="GGtools"))
pp[[1]][1:3]
pp[[2]]
personalHap(pp)
}
\keyword{ models }
