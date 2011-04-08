\name{relocate}
\alias{relocate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assist in the transport between systems of ff data managed by GGtools 
}
\description{
assist in the transport between systems of ff data managed by GGtools 
}
\usage{
relocate(old, new, obj, ffind = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{old}{
string to be replaced in the physical filename attribute on old system
}
  \item{new}{
string to be substituted for \code{old} in the physical filename attribute on old system
}
  \item{obj}{
manager object
}
  \item{ffind}{
index of file in fflist to be altered
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a new manager instance
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (old, new, obj, ffind = 1) 
{
    inst = obj@fflist[[ffind]]
    fref = attr(attributes(inst)[["physical"]], "filename")
    ans = gsub(old, new, fref)
    attr(attributes(inst)[["physical"]], "filename") = ans
    obj@fflist[[ffind]] = inst
    obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
