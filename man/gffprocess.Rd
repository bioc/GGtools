\name{gffprocess}
\alias{gffprocess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transform a collection of gff3 into a single tabix-indexed gff3
}
\description{
process a collection of gff3 into a single tabix-indexed gff3 using
unix utilities to minimize memory requirements
}
\usage{
gffprocess(basename = "fullyri100k", n_in = 44, headpatt = "_1A", tmpForSort = "/freshdata/tmp")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{basename}{
basename of the resulting .gff3.gz(.tbi) file
}
  \item{n_in}{
number of gff3 files to be processed -- used for consistency check
against length(dir(patt="gff3")) 
}
  \item{headpatt}{
pattern to identify file for the 'top' gff3 to be used as the
contents are concatenated
}
  \item{tmpForSort}{
name of a folder that unix sort will use as a temporary directory
}
}
\details{
The purpose of this utility is to exploit unix shell tools
to unify a collection of gff3 files generated using \code{link{All.cis}}.
The use case is cluster-based per-chromosome (or split chromosome)
cis-testing generating a large number of GRanges that are transformed
to gff3 to allow targeted interrogation.
}
\value{
Used for side effects.  Will fail if any unix utility call
via system() returns nonzero value.  Returns NULL otherwise.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
