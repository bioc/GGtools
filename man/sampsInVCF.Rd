\name{sampsInVCF}
\alias{sampsInVCF}
\alias{TabixFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
enumerate samples available in a VCF file
}
\description{
enumerate samples available in a VCF file
}
\usage{
sampsInVCF(tf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tf}{
instance of \code{\link[Rsamtools]{TabixFile}} referring
to a tabix-indexed VCF
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
\note{
This package exports \code{\link[Rsamtools]{TabixFile}} for the sake of the example below.
}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
tf = TabixFile(system.file("vcf/CEU.exon.2010_09.genotypes.vcf.gz", package="GGtools"))
sampsInVCF(tf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
