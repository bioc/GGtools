\name{pifdr}
\alias{pifdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
utility for computing plug-in FDR
}
\description{
utility for computing plug-in FDR
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pifdr( obs, perms, npts=1999, applier=sapply)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
observed association scores
%%     ~~Describe \code{obs} here~~
}
  \item{perms}{
vector of association scores under permutation; length should
be integer multiple of \code{length(obs)}
%%     ~~Describe \code{ps} here~~
}
  \item{npts}{number of points spanning the range
of \code{obs} to be used for a lossy grid-based computation;
only used of \code{length(obs)>npts}.
  \item{applier}{
function that iterates the computation
%%     ~~Describe \code{applier} here~~
}
}
\details{
As currently implemented the algorithm is quadratic in \code{length(obs)}.
While it is possible to get a unique FDR value for every element of \code{obs},
an approximate approach yields practically identical precision and by default this
will be used for \code{obs} with length 2000 or more.  In this case,
\code{\link{approx}} is used with \code{rule=2} to interpolate from the
grid-based FDR estimates back to the data values.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
vector of plug-in FDR estimates congruent to \code{obs}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Hastie Tibshirani and Friedman Elements of Statistical
Learning ch 18.7
%% ~put references to the literature/web site here ~
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
X = rchisq(100,1)
Y = rchisq(300,1)
qqplot(pifdr(X,Y), rchisq(100,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
