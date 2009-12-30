\name{snpm2phase}
\alias{snpm2phase}
%\alias{parsePhPairs}
\alias{parsePh.out}
\alias{personalHap}
%\alias{invokePhase}
\alias{phaseInput-class}
\alias{show,phaseInput-method}
%\alias{invokePhase,phaseInput,chrnum,character,character,character,logical-method}
%\alias{invokePhase,snp.matrix,chrnum,character,character,character,logical-method}

\title{convert information in a snp.matrix to PHASE input format;
invokePhase can run a suitably installed version of PHASE  }
\description{convert information in a snp.matrix to PHASE input format;
invokePhase can run a suitably installed version of PHASE  }

\usage{
snpm2phase(snpm, cnum, outfilename)
parsePh.out(fn)
personalHap(x)
}
\arguments{
  \item{snpm}{ snp.matrix instance}
  \item{cnum}{ chromosome number as chrnum instance}
  \item{outfilename}{ character name of file to write}
  \item{fn}{ character name of PHASE .out file to read}
  \item{x}{ output of parsePh.out}
}
\details{ follows phase 2.1 documentation for input format

a phaseInput container class can store relevant metadata
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

class(pin)
getClass("phaseInput")
pin
readLines(tt)
unlink(tt)
pp = parsePh.out(system.file("phaseOut/cpne1_20k.out", package="GGtools"))
pp[[1]][1:3]
personalHap(pp)
}
\keyword{ models }
