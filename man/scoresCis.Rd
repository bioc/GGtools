\name{scoresCis}
\alias{scoresCis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visualize a gene model with cis-eQTL association scores (-log FDR by default) on the
basis of a ciseqByCluster data.table output
}
\description{
visualize a gene model with cis-eQTL association scores (-log FDR by default) on the
basis of a ciseqByCluster data.table output
}
\usage{
scoresCis(sym = "ORMDL3", cisRun, 
   cisannopk = "lumiHumanAll.db", radius = 1e+05, pad = 1000, 
   txScore = function(x) -log10(x + (1e-05)), ylim = c(0, 4), 
   genometag = "hg19", plot.it = TRUE, laxistag = "-log10 FDR: ", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sym}{
gene symbol to be resolved into probe id using cisannopk
}
  \item{cisRun}{
data.table output of ciseqByCluster 
}
  \item{cisannopk}{
Annotation resource, often a ChipDb instance
}
  \item{radius}{
radius to be added to gene model for display, should typically agree with that used in the search
}
  \item{pad}{
some extra space
}
  \item{txScore}{
function that will transform fdr for rendering
}
  \item{ylim}{
vertical limits for fdr display
}
  \item{genometag}{
coordinates from this build of genome
}
  \item{plot.it}{
logical dictating whether plotTracks will be run
}
  \item{laxistag}{
token used to tell what units are used on vertical axis
}
  \item{\dots}{
not used
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a list of Gviz tracks, invisibly returned
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
%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The Bioconductor workflow on cloud-enabled cis-eQTL analysis.
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
