\name{chr20meta}
%c15GGceu.rda    c7GGceu.rda     chr7meta.rda    racHMceu20.rda
%c20GG.rda       chr15meta.rda   geneLocs.rda    racHMceu7.rda
%c20GGceu.rda    chr20meta.rda   racHMceu15.rda
%
\alias{chr20meta}
\alias{chr7meta}
\alias{chr15meta}
\alias{c20gt}
\title{representations of HapMap snp metadata}
\description{
representations of HapMap snp metadata
}
\usage{
data(chr20meta)
data(c20gt)
}
\details{
The ``meta'' structures are instances of class \code{\link{data.frame}},
created from initial columns of HapMap chromosome-specific
genotyping results for CEPH CEU participants.  The
\code{c20gt} is a named matrix of unphased genotyping results.

}
\value{
instances of class \code{\link{data.frame}}
}
\references{ 
HapMap project
 }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(chr20meta)
# be patient
chr20meta[1:4,]
data(c20gt)
dim(c20gt)
c20gt[1:4,1:4]
}
\keyword{ models }
