\name{ggreg}
\alias{ggreg}
\title{regression analysis for
additive effect of genotype on expression  }
\description{regression analysis for
additive effect of genotype on expression
}
\usage{
ggreg(gge, psid, snpinds=NULL, model=c("additive", "general")[1])
}
\arguments{
  \item{gge}{ \code{\link{ggExprSet}} instance}
  \item{psid}{ name of the dependent
variable in the regression analysis,
typically identified using an affy probe set id (more generally,
rownames element for the exprs response for
the gge argument }
  \item{snpinds}{ numerical indices of snps to
be analyzed; if NULL, all will be.}
  \item{model}{character tag indicating selection of additive model
(number of rare alleles regarded as a continuous predictor) or
general model (separate effects estimated for homozygous 
rare and heterozygous)}
}
\value{a list of \code{\link{lm}} objects
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(c20GGceu)
ggcpne1 = ggreg(c20GGceu, "206918_s_at", which(names(snps(c20GGceu)) ==
  "rs6060535"))
summary(ggcpne1[[1]])
}
\keyword{ models }
