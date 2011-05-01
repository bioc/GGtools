\name{mrefhap2sm}
\alias{mrefhap2sm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transform MACH-supplied haplotype data for imputation into a SnpMatrix instance
}
\description{
transform MACH-supplied haplotype data for imputation into a SnpMatrix instance
}
\usage{
mrefhap2sm(gzfn, snpids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gzfn}{
name of gzipped file with haplotype sequences
}
  \item{snpids}{
vector of unique SNP ids for the haplotype elements
}
}
\details{
uses \code{\link[snpStats]{read.snps.long}}.  The MACH group provides
haplotypes as two long strings of nucleotide codes per individual.
}
\value{
an instance of \code{\link[snpStats]{SnpMatrix-class}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
smhapf = system.file("machHap/c20small.hap.gz", package="GGtools")
snidf = gzfile(system.file("machHap/chr20.snps.gz", package="GGtools"))
snids = scan(snidf, "")
sm = mrefhap2sm( smhapf, snids )
sm
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
