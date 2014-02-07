\name{qqhex}
\alias{qqhex}
\alias{binqq}
\alias{binnedQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
obtain qqplot coordinates for
the specific case of comparing a given distribution
to that of multiple realizations from
a permutation distribution, and bin these coordinates using hexbin, useful
for very large samples
}
\description{
obtain qqplot coordinates for
the specific case of comparing a given distribution
to that of multiple realizations from
a permutation distribution, and bin these coordinates using hexbin, useful
for very large samples
}
\usage{
binnedQQ(dt, nxbins=20,
  ylim=c(0,76), xlim=c(0,30), end45=5, thrs=c(0,.001,.005,.01,.05), 
  tempmar = c(6,4,4,5), ...)
qqhex(sco, p1, p2, p3, fdr, nxbins = 20, thrs = c(0, 0.001, 0.005, 0.01, 0.05))
binqq(qqob, ylim = c(0, 76), xlim = c(0, 30), end45=5, ...)
}
\arguments{
  \item{dt}{a \code{\link[data.table]{data.table}} instance
with association scores and scores obtained under permutation along
with FDR, as returned by \code{\link{cgff2dt}} or \code{\link{ciseqByCluster}}}
  \item{sco}{
numeric vector of observed statistics
}
  \item{p1}{
realization of null distribution for \code{sco}, independent of \code{p2} and
\code{p3}
}
  \item{p2}{
realization of null distribution for \code{sco}, independent of \code{p1} and
\code{p3}
}
  \item{p3}{
realization of null distribution for \code{sco}, independent of \code{p1} and
\code{p2}
}
  \item{fdr}{
vector of FDR associated with elements of \code{sco}
}
  \item{nxbins}{
number of bins to be used for samples from the null distribution
}
  \item{thrs}{
vector of thresholds in FDR to be used for ruling the plot
}
  \item{qqob}{
for binhex(), output of qqhex
}
  \item{ylim}{
vertical limits of rendering
}
  \item{xlim}{
horizontal limits of rendering
}
  \item{end45}{
a segment is drawn from (0,0) to (end45,end45) to depict the line of identity
}
  \item{tempmar}{
numerical vector with 4 elements serving as a temporary setting of
the \code{mar} parameter of \code{\link[graphics]{par}}
}
  \item{\dots}{
not currently used
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
for qqhex, a list with elements
\item{hb }{output of \code{\link{hexbin}}}
\item{thrs }{vector of input \code{thrs}}
\item{scothrs }{vector of observed statistics corresponding to FDRs in \code{thrs}}
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
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
opar = par(no.readonly=TRUE)
set.seed(123)
x = c(rchisq(9000,1), rchisq(1000,12))
nn = lapply(1:3, function(x) rchisq(10000,1))
fd = pifdr(x, unlist(nn))
qqh = qqhex(x, nn[[1]], nn[[2]], nn[[3]], fd)
par(mar=c(4,4,4,7))
binqq(qqh,xlim=c(0,10), ylim=c(0,20))
mtext(4, "FDR")
par(opar)
}
\keyword{ models }
