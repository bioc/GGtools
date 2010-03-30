\name{snp130locs}
\alias{snp130locs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
prototypical function for creation of IRanges-based SNP location data 
}
\description{
prototypical function for creation of IRanges-based SNP location data 
}
\usage{
snp130locs(chr, start, end)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
scalar string with prefix "chr"
}
  \item{start}{
numeric start value, typically 1
}
  \item{end}{
numeric end value, typically length in bases of chromosome
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a \code{\link[rtracklayer]{ucscTableQuery}} output
%%%  ~Describe the value returned
%%%  If it is a LIST, use
%%%  \item{comp1 }{Description of 'comp1'}
%%%  \item{comp2 }{Description of 'comp2'}
%%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## The function is currently defined as
function (chr, start, end) 
{
    sess = browserSession()
    quer = ucscTableQuery(sess, "snp130", GenomicRanges(start, 
        end, chr))
    tableName(quer) = "snp130"
    track(quer)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
