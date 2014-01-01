\name{collectBest}
\alias{collectBest}
\alias{collectFiltered}
\alias{cis.FDR.filter.best}
\alias{cis.FDR.filter.SNPcentric}
\alias{cis.FDR.filter.SNPcentric.complete}
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
 targetname = "harvest", 
 mafs = c(0.01, 0.02, 0.025, 0.03333, 0.05, 0.075, 0.1), 
 hidists = c(10000, 25000, 50000, 75000, 1e+05, 250000), interimSaves=FALSE)

collectFiltered( fns, targetname="harvest",
   mafs = c(.01, .02, .025, .03333, .05, .075, .1),
   hidists = c(10000, 25000, 50000, 75000, 100000, 250000),
   filterFun = cis.FDR.filter.best, filtApplier=lapply,
   interimSaves=FALSE) 

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
  \item{filterFun}{function like GGtools:::cis.FDR.filter.best}
  \item{filtApplier}{function like lapply}
  \item{interimSaves}{logical, if TRUE save list at each maf/dist transition}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
#
# contents of fns are two chromsomes of cis runs for CEU
#
  fns = dir(system.file("rdas", package="GGtools"), full=TRUE)
  cc = collectBest(fns, mafs=c(.01, .05), hidists=c(10000, 50000))
  sapply(cc, sapply, function(x) sum(x$fdr <= 0.01))
#
# this tells us which to keep
#
  kp = cc[["0.05"]][["50000"]]
  kp = kp[kp$fdr <= 0.01,]
#
# the hits are in the table above; the following function
# retrieves the initial scores giving rise to the filtered
# hits
#
  pullHits = function(fns, atts) {
    tmp = lapply(fns, function(x) get(load(x)))
    kl = lapply(tmp, function(x) paste(names(x), x$snp, sep=":"))
    attk = paste(atts$genes, atts$bestsnp, sep=":")
    tmp = lapply(1:length(tmp), function(x) tmp[[x]][ match( attk, kl[[x]], nomatch=0 ) ])
    curans = do.call(c, lapply(tmp, as, "GRanges"))
    neword = match( attk, paste(names(curans), curans$snp, sep=":"))
    newfdr = atts$fdr[neword]
    curans$fdr = newfdr
    curans
    }
  pullHits( fns, kp )
#
#
#
#   after executing code in example for All.cis (protected by dontrun)
#   and running save(f1, file="f1.rda"), the following will work
# genewise max score
    cf1 = collectFiltered("f1.rda", mafs=.02, hidists=25000, targetname="gwise")
# SNPwise scores, all
    cf2 = collectFiltered("f1.rda", mafs=.02, hidists=25000, targetname="swise",
         filterFun = cis.FDR.filter.SNPcentric.complete )
# SNPwise scores, best per SNP when SNP is cis to multiple genes
    cf3 = collectFiltered("f1.rda", mafs=.02, hidists=25000, targetname="swise2",
         filterFun = cis.FDR.filter.SNPcentric )
  }  # end dontrun

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
