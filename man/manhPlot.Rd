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
manhPlot(probeid, mgr, ffind, namedlocvec = NULL, locGRanges = NULL, plotter = smoothScatter, tx = function(x) -log10(1 - pchisq(x, 1)), xlab = paste("pos. on ", names(fflist(mgr))[ffind]), ylab = "-log10 p", ...)
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
\dontrun{
if (!exists("e1")) example(eqtlTests)  # creates e1, cptag, c20r
# use ffind=1 below because you have confined attention to chr20
manhPlot( cptag, e1, ffind=1, locGRanges=c20r, cex=3)
setwd(curd)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
