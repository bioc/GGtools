\name{calfig}
\alias{calfig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot results of appraise to exhibit model calibration
}
\description{
plot results of appraise to exhibit model calibration
}
\usage{
calfig(colist, 
 tabs, ind = 10, hfudgetxt = 0.0155, tickend = 0.16, 
 tickgap = 0.02, ylimin = c(-0.01, 0.16), 
 xlimin = c(-0.01, 0.16), fraccex = 0.8, fuselast = 0)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{colist}{\code{appraise} output object loaded
from *coeflist.rda, a list of \code{summary(biglm.obj)$mat} for
different models}
  \item{tabs}{\code{appraise} output object loaded
from *tabs.rda, a list of tables counting SNPs in bins
of predicted probabilities}
  \item{ind}{index or name of model to be plotted}
  \item{hfudgetxt}{distance to move rendered fractions relative
to bin x coordinate}
  \item{tickend}{
maximum value at which axis tick mark will be plotted
}
  \item{tickgap}{
axis will have ticks at \code{seq(0, tickend, tickgap)}
}
  \item{ylimin}{
ylim setting for rendering
}
  \item{xlimin}{
xlim setting for rendering
}
  \item{fraccex}{
cex setting for fraction rendering
}
  \item{fuselast}{
if data are sparse in entries of high predicted probability,
you can fuse the nearby cells up to the end -- pick fuselast=2
for final 2 cells, 3 for final 3 and so on
}
}
\details{
can be finicky ... assumes \code{\link{appraise}} used in a fairly vanilla way
}
\value{
renders
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
