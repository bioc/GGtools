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
builds a \code{\link{racExSet-class}} instance
}
\value{
builds a \code{\link{racExSet-class}} instance
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(chr20GGdem)
fi = system.file("fileDemos/geno20demo.txt.gz", package="GGtools")
HMworkflow(fi, exprs(chr20GGdem)[1:100,],
  phenoData(chr20GGdem), experimentData(chr20GGdem),
  "hgfocus")
}
\keyword{ models }
