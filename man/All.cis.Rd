\name{All.cis}
\alias{All.cis}
\alias{chrFilter}
\alias{mcwAllCis-class}
\alias{show,mcwAllCis-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function that computes score tests for all SNP cis to genes, with
flexible filtering
}
\description{
function that computes score tests for all SNP cis to genes, with
flexible filtering
}
\usage{
All.cis( config = new("CisConfig"), ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{config}{
instance of class \code{\link{CisConfig-class}}
}
  \item{\dots}{
passed to eqtlTests
}
}
\details{
returns score statistics for assocations of all SNP cis to genes,
in a GRanges instance, with range names given by probes; metadata supplied
SNP location, name, and score
}
\value{
instance of \code{\link{cisRun-class}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
  cc = new("CisConfig")
  chrnames(cc) = "21"
  lkp = try(library(parallel))
  if (!inherits(lkp, "try-error")) {
     nc = min(10, detectCores())
     options(mc.cores=nc)
     geneApply(cc) = mclapply
     }
  estimates(cc) = FALSE
  unix.time(f1 <- All.cis( cc ))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
