\name{sampsInVCF}
\alias{sampsInVCF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
enumerate samples available in a VCF file
}
\description{
enumerate samples available in a VCF file
}
\usage{
sampsInVCF(tf, chrToTry = "1", poke = IRanges(1, 1000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tf}{
instance of \code{\link[Rsamtools]{TabixFile}} referring
to a tabix-indexed VCF
}
  \item{chrToTry}{
character label for chromosome
}
  \item{poke}{
region to query
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
vector of available sample identifiers
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
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tf, chrToTry = "1", poke = IRanges(1, 1000)) 
{
    v1 = ScanVcfParam(fixed = "ALT", info = NA, geno = "GT", 
        which = GRanges("1", IRanges(1, 1000)))
    pro = readVcf(tf, genome = "hg19", param = v1)
    colnames(geno(pro)[[1]])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
