\name{All.cis}
%\alias{cisAssoc}
\alias{cisScores}
\alias{All.cis}
\alias{chrFilter}
\alias{mcwAllCis-class}
\alias{show,mcwAllCis-method}
\alias{add878}
\alias{addgwhit}
\alias{inflammFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
functions that compute score tests for all SNP cis to genes, with
flexible filtering
}
\description{
function that computes score tests for all SNP cis to genes, with
flexible filtering
}
\usage{
cisScores( config = new("CisConfig"), ... )
All.cis( config = new("CisConfig"), ... )
addgwhit(ans, traitFilter=force, vname="isgwashit") 
add878(ans)
inflammFilter(gwtagger)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{config}{
instance of class \code{\link{CisConfig-class}}
}
  \item{\dots}{
passed to eqtlTests
}
%  \item{summex}{
%instance of \code{\link[GenomicRanges]{SummarizedExperiment-class}}
%}
%  \item{vcf.tf}{
%instance of \code{\link[Rsamtools:TabixFile-class]{TabixFile}}, pointing to a
%tabix-indexed VCF file
%}
%  \item{rhs}{
%headless formula instance; covariates can be identified
%in colData(summex) 
%}
%  \item{nperm}{
%number of permutations of assay against genotype to be retained,
%should be a small number, say 3, unless the feature set is very
%small
%}
%  \item{cisradius}{
%number of bp up- and down-stream  of each rowRanges(summex) range
%to be searched for variants in \code{vcf.tf}
%}
%  \item{stx}{
%not used, intended for future coding as a filter for assay data
%}
%  \item{vtx}{
%not used, intended for future coding as a filter for VCF data
%}
%  \item{snfilt}{
%seqnames filter to harmonize differences between chromosome names
%in summex and vcf.tf -- defaults to 'remove chr prefix'
%}
%  \item{genome}{
%conventional tag like 'hg19'
%}
%  \item{assayind}{
%numeric, take quantitations from \code{assays(summex)[[assayind]]}
%}
%  \item{lbmaf}{
%lower bound on MAF of variants to be retained
% }
  \item{ans}{cisRun-like entity to which additional annotation
    will be bound by addgwhit or add878}
  \item{gwtagger}{GRanges like gwastagger in gwascat data elements}
  \item{traitFilter}{function that returns a gwastagger-like GRanges, see
\code{inflammFilter}}
  \item{vname}{name to be used for new data.table column added by
addgwhit}
}
\details{
cisScores (All.cis) 
returns score statistics for assocations of all SNP cis to genes,
in a GRanges instance, with range names given by probes; metadata supplied
SNP location, name, and score

cisAssoc targets SummarizedExperiment instances for molecular
phenotype measures and VCF for variant data

\code{addgwhit} and \code{add878} will use
GWAS hit information or ChromHMM labeling to annotation
ranges
}
\value{
for cisScores: instance of \code{\link{cisRun-class}}

for cisAssoc: a GRanges with information on observed and
permuted test scores per locus/feature pair
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
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
   system.time(f1 <- cisScores( cc ))
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
 system.time(f2 <- cisScores( cc ))
#
#
 inflammFilter # to make more restrictive predicate for prediction
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
