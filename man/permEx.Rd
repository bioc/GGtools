\name{permEx}
\alias{permEx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
permute expression data against genotype data in an smlSet
}
\description{
permute expression data against genotype data in an smlSet
}
\usage{
permEx(sms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sms}{
an instance of \code{\link{smlSet-class}}
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
an instance of \code{\link{smlSet-class}}
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (!exists("hmceuB36.2021")) data(hmceuB36.2021)
library(illuminaHumanv1.db)
cptag = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
indc = which(featureNames(hmceuB36.2021) == cptag[1])
hm = hmceuB36.2021[c(indc,1:19),]  # reduce problem
td = tempdir()
curd = getwd()
setwd(td)
time.lapply = unix.time(e1 <- eqtlTests( hm, ~male ))
e1
hmp = permEx(hm)
e1perm = eqtlTests(hmp, ~male, targdir="permfoo", runname="permrun")
topFeats(probeid=cptag, mgr=e1, ffind=1, anno="illuminaHumanv1.db", useSym=FALSE)
topFeats(probeid=cptag, mgr=e1perm, ffind=1, anno="illuminaHumanv1.db", useSym=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
