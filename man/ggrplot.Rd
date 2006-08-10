\name{ggrplot}
\alias{ggrplot}
\title{plot a regression model for
additive effect of genotype on expression  }
\description{plot a regression model for
additive effect of genotype on expression
}
\usage{
ggrplot(gge, genesym, snpname=NULL, plat="hgfocus", psid=NULL)
}
\arguments{
  \item{gge}{ instance of \code{\link{ggExprSet}}}
  \item{genesym}{ HUGO gene symbol}
  \item{snpname}{ dbSNP rs number for SNP}
  \item{plat}{ platform name for symbol resolution}
  \item{psid}{ probe set id in case there are multiple
probe sets for the given HUGO symbol}
}
%\details{
%}
\value{invisibly returns list with dependent and independent
variables and the \code{\link{lm}} object representing
the regression
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(c20GGceu)
library(hgfocus)
ggrplot(c20GGceu, "CPNE1", "rs6060535")
}
\keyword{ models }
