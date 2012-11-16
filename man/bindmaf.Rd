\name{bindmaf}
\alias{bindmaf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bind testing metadata to a best.cis.eQTLs result
}
\description{
bind testing metadata to a best.cis.eQTLs result

}
\usage{
bindmaf(smpack = "GGdata", smchr = "20", obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smpack}{
name of a package to which \code{\link[GGBase]{getSS}} can be applied
to generate an instance of \code{\link[GGBase]{smlSet-class}}
}
  \item{smchr}{
the chromosome name as used in the names of the \code{smList} output for
the \code{getSS} result
}
  \item{obj}{
an instance of \code{\link{mcwBestCis-class}} generated using the
package named in \code{smpack}
}
}
\details{
computes the MAF of most highly associated SNP per gene, and distance between
that SNP and the transcription limits of the gene, assigning 0 for this if the
SNP lies within the transcription limits
}
\value{
a GRanges instance
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
\note{
This will be used to stratify the permuted scores.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
 b1 = best.cis.eQTLs(chr="20")  # sharply filtered
 b1b = bindmaf(obj=b1)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
