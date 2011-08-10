\name{manhPlot}
\alias{manhPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
manhattan plot for an eqtlTests result
}
\description{
manhattan plot for an eqtlTests result
}
\usage{
manhPlot(probeid, mgr, ffind, namedlocvec = NULL, locGRanges = NULL, plotter = smoothScatter, tx = function(x) -log10(1 - pchisq(x, 1)), xlab = paste("pos. on ", names(fflist(mgr))[ffind]), ylab = "-log10 p", suppressGeneLoc=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probeid}{
element of colnames of fflist(mgr)[[ffind]] -- the gene of interest, typically
}
  \item{mgr}{
an instance of eqtlTestsManager
}
  \item{ffind}{
index of the ff file of interest -- typically identifying a chromosome
where SNP locations define the x-axis of the plot
}
  \item{namedlocvec}{
a vector with named elements, giving SNP locations
}
  \item{locGRanges}{
a GRanges instance with SNP locations
}
  \item{plotter}{
function to be used for rendering
}
  \item{tx}{
the numbers acquired from the manager are assumed to be chi-squared(1) -- this
function can be altered to define how the y axis is derived from manager contents
}
  \item{xlab}{
label for x axis
}
  \item{ylab}{
label for y axis
}
  \item{suppressGeneLoc}{
logical; if true, will refrain from trying to indicate gene location
on plot.  Important to have TRUE when a trans association is being plotted.
}
  \item{\dots}{
passed to plotting function
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%%  ~Describe the value returned
%%%  If it is a LIST, use
%%%  \item{comp1 }{Description of 'comp1'}
%%%  \item{comp2 }{Description of 'comp2'}
%%% ...
%}
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
if (require(SNPlocs.Hsapiens.dbSNP.20100427)) {
 library(GGdata)
 hm2ceuSMS = getSS("GGdata", "20", renameChrs="chr20")
 library(illuminaHumanv1.db)
 cptag = get("CPNE1", revmap(illuminaHumanv1SYMBOL))
 indc = which(featureNames(hm2ceuSMS) == cptag[1])
 hm = hm2ceuSMS[c(indc,1:19),]  # reduce problem
 hm = hm[chrnum("chr20"),]
 td = tempdir()
 curd = getwd()
 setwd(td)
 e1 <- eqtlTests( hm, ~male, targdir="mplex" )
 c20  = getSNPlocs("ch20", as.GRanges=TRUE)
 sr = ranges(c20)
 sr = GRanges(seqnames="chr20", sr)
 elementMetadata(sr) = elementMetadata(c20)
 names(sr) = paste("rs", elementMetadata(sr)$RefSNP_id, sep="")
# use ffind=1 below because you have confined attention to chr20
 manhPlot( cptag, e1, ffind=1, locGRanges=sr, cex=3)
 setwd(curd)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
