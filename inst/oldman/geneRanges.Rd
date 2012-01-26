\name{geneRanges}
\alias{geneRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute a named GRanges instance with extents of a set of genes
on a chromosome
}
\description{
compute a named GRanges instance with extents of a set of genes
on a chromosome
}
\usage{
geneRanges(genomeOrPkgOrTxDb, chr = NULL, is.annopkg = FALSE, 
  extend = 0, applier=lapply)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genomeOrPkgOrTxDb}{
A string or TranscriptDb instance.  if \code{is.annopkg} is \code{TRUE}
a string must be passed and must name an installed chip or org. annotation package. 
If a string and \code{is.annopkg} is \code{FALSE}, must be either
\code{"hg18"} or \code{"hg19"} in which case a TxDb.Hsapiens package
will be used.  Otherwise a \code{TranscriptDb} instance must be passed.
}
  \item{chr}{
a string naming the chromosome for which extents are requested.
}
  \item{is.annopkg}{
logical taking value \code{TRUE} to indicate that \code{genomeOrPkgOrTxDb}
names an annotation package.
}
  \item{extend}{
numeric indicating radius of symmetric extension desired around annotated extent
}
  \item{applier}{
application function used to obtain extents of genes when iterating over transcriptsBy
result -- this seemed slow for an ENSEMBL TxDb so multicore application can be considered
by setting this parameter to mclapply.
}
}
\details{
This function is intended to help users get locations for named genes.
The CHRLOC mappings are used for org. or chip annotation packages.
The degapped extents of transcriptsBy results are returned otherwise.
}
\value{
A names GRanges instance.
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
This function will not be exported in GGtools 4.x.
}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
if (require(org.Hs.eg.db)) {
  geneRanges("org.Hs.eg.db", is.annopkg=TRUE, chr="chr20")
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
