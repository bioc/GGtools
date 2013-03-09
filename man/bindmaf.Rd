\name{bindmaf}
\alias{bindmaf}
\alias{meta.bindmaf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bind testing metadata to a best.cis.eQTLs result
}
\description{
bind testing metadata to a best.cis.eQTLs result

}
\usage{
bindmaf(smpack = "GGdata", smchr = "20", obj, SSgen=GGBase::getSS)
meta.bindmaf (smpackvec=c("GGdata", "hmyriB36"),
     smchr="20", obj, usemaxMAF=FALSE, SSgen=GGBase::getSS) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smpack}{
name of a package to which \code{\link[GGBase]{getSS}} can be applied
to generate an instance of \code{\link[GGBase]{smlSet-class}}
}
  \item{smpackvec}{a vector of candidate smpack values for metaanalysis
across populations or tissues}
  \item{smchr}{
the chromosome name as used in the names of the \code{smList} output for
the \code{getSS} result
}
  \item{obj}{
an instance of \code{\link{mcwBestCis-class}} generated using the
package named in \code{smpack}
}
  \item{usemaxMAF}{if TRUE, label a SNP with maximum MAF observed across
populations, otherwise compute the MAF for the combined genotypes across
populations represented by the various smlSet instances generated with
the \code{smpackvec} spec.
  }
 \item{SSgen}{function to be used to create smlSet instance
for testing -- in general, GGBase::getSS has been used to pull
the ExpressionSet and SnpMatrix data from a named package, but
in some cases a specialize task is needed to create the desired
smlSet.  Whatever is passed to SSgen must return an smlSet instance.}
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
