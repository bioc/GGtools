\name{All.cis}
\alias{cisScores}
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
cisScores( config = new("CisConfig"), ... )
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
   genome(cc) = "hg19"
   lkp = try(library(parallel))
   if (!inherits(lkp, "try-error")) {
      nc = min(10, detectCores())
      options(mc.cores=nc)
      geneApply(cc) = mclapply
      }
   estimates(cc) = FALSE
   set.seed(1234)
   unix.time(f1 <- cisScores( cc ))
 #
 # demonstrate adding annotation on chromatin state and gwas status
 #
 eprops = function(ans) {
 #
 # only adds fields to values() of the input
 #
  data(hmm878)
  ac = as.character
  eqr = GRanges(ac(seqnames(ans)), IRanges(ans$snplocs, width=1))
  fo = findOverlaps(eqr, hmm878)
  chromcat878 = factor(rep("none", length(ans)), levels=c(unique(hmm878$name), "none"))
  chromcat878[ queryHits(fo) ] = factor(hmm878$name[subjectHits(fo)])
  ans$chromcat878 = chromcat878
 
  if (require(gwascat)) {
    data(gwastagger)
    isgwashit = 1*(overlapsAny(eqr, gwastagger) | ans$snp %in% gwastagger$tagid) # allow match by loc or name
    ans$isgwashit = isgwashit
    }
  ans
 }
 extraProps(cc) = eprops
 set.seed(1234)
 unix.time(f2 <- cisScores( cc ))
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
