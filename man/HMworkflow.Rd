\name{HMworkflow}
\alias{HMworkflow}
\title{function to bind together
HapMap genotyping results and expression data  }
\description{ 
function to bind together
genotyping results and expression data
}
\usage{
HMworkflow(gzfn, emat, pd, mi, anno)
}
\arguments{
  \item{gzfn}{ path to gzipped genotyping results }
  \item{emat}{ matrix of conforming expression results
(colnames match some of those in gzfn) }
  \item{pd}{ \code{\link[Biobase]{phenoData}} instance }
  \item{mi}{ \code{\link[annotate]{MIAME}} instance}
  \item{anno}{ string indicating platform, such as
\code{"hgfocus"}}
}
\details{
builds a \code{\link{racExSet}} instance
}
\value{
builds a \code{\link{racExSet}} instance
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(chr20GGceu)
fi = system.file("fileDemos/geno20demo.txt.gz", package="GGtools")
HMworkflow(fi, exprs(chr20GGceu)[1:100,],
  phenoData(chr20GGceu), experimentData(chr20GGceu),
  "hgfocus")
}
\keyword{ models }
