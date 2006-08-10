\name{regseq}
\alias{regseq}
\title{ compute a sequence of regression models
relating genotype to expression  }
\description{
compute a sequence of regression models
relating genotype to expression  
}
\usage{
regseq(gge, psid, snpinds, snpmeta, ...)
}
\arguments{
  \item{gge}{ instance of \code{\link{ggExprSet}}}
  \item{psid}{ probe set id (geneNames element)}
  \item{snpinds}{ numerical indices of snps to be used (relative
to order in the snps component of \code{gge}}
  \item{snpmeta}{ a data.frame with position information
(see \code{\link{chr20meta}} for example)}
  \item{\dots}{additional arguments (such as \code{model="general"}
    to be passed to \code{\link{ggreg}}.}
}
\details{
}
\value{
}
\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{  }


\seealso{  }

\examples{
# used by detScreen
}
\keyword{ models }
