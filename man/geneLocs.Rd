\name{geneLocs}
%c15GGceu.rda    c7GGceu.rda     chr7meta.rda    racHMceu20.rda
%c20GG.rda       chr15meta.rda   geneLocs.rda    racHMceu7.rda
%c20GGceu.rda    chr20meta.rda   racHMceu15.rda
%
\alias{geneLocs}
\title{gene metadata from NCBI}
\description{
gene metadata from NCBI
}
\usage{
data(geneLocs)
}
\details{
A \code{\link{data.frame}} instance with information relating
gene symbol to chromosome and providing start and end positional
information.
}
\value{
instance of class \code{\link{data.frame}}
}
\references{ 
NCBI
 }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(geneLocs)
# be patient
geneLocs[1:4,]
}
\keyword{ models }
