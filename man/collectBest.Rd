\name{collectBest}
\alias{collectBest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
given a collection of All.cis outputs (cisRun instances)
compute FDRs for various filterings
}
\description{
given a collection of All.cis outputs (cisRun instances)
compute FDRs for various filterings
}
\usage{
collectBest(fns, 
 targetname = "colls", 
 mafs = c(0.01, 0.02, 0.025, 0.03333, 0.05, 0.075, 0.1), 
 hidists = c(10000, 25000, 50000, 75000, 1e+05, 250000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fns}{
names of .rda with the cisRun outputs
}
  \item{targetname}{
basename of rda file to be emitted
}
  \item{mafs}{
lower bounds on MAF for filtering
}
  \item{hidists}{
upper bounds on cis radius for filtering
}
}
\details{
\code{\link{pifdr}} is repeatedly used to generate
conditional plugin FDR for different filtering criteria
}
\value{
a list of lists is written to disk incrementally, as the job
can be long running
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
This is the workhorse of sensitivity analysis.  Permits counting
of genes with eQTL at selected FDR for various criteria on cis
radius and MAF.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
