\name{dropMonomorphies}
\alias{dropMonomorphies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
remove monomorphic loci from an smlSet instance
}
\description{
remove monomorphic loci from an smlSet instance
}
\usage{
dropMonomorphies(sms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sms}{
instance of \code{\link[GGBase]{smlSet-class}}
}
}
\details{
uses \code{\link[snpStats]{col.summary}} to determine
monomorphy
}
\value{
instance of \code{\link[GGBase]{smlSet-class}}
}
\examples{
library(GGdata)
gg20 = get("GGdata", "20")
dim(smList(gg20)[[1]])
dim(smList(dropMonomorphies(gg20))[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
