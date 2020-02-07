\name{getCisMap}
\alias{getCisMap}
\alias{show,cisMap-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create, using
Bioconductor annotation resources, a structure that enumerates SNP in the vicinity of ('cis' to) genes
}
\description{
create a structure that enumerates SNP in the vicinity of ('cis' to) genes
}
\usage{
getCisMap(radius = 50000, gchr = "20", 
  schr = "ch20", geneannopk = "illuminaHumanv1.db", 
  snpannopk = snplocsDefault(),
  as.GRangesList = FALSE, excludeRadius=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{radius}{
How far, in bases, up or down stream from the asserted coding region limits
to include SNP
}
  \item{gchr}{
the token to be used to acquire locations for probes on a
specified chromosome, using
\code{revmap([dbpk]CHR)}
}
  \item{schr}{
the token to be used to acquire locations for SNP on a specified
chromosome, using \code{\link[BSgenome:SNPlocs-class]{snplocs}}
}
  \item{geneannopk}{
character string naming a Bioconductor .db expression chip annotation package;
or a complex string with first part naming a Bioconductor FDb.* annotation
package, colon separator, and a second string naming the getter hook that
when called returns a GRanges with names corresponding to features and ranges
corresponding to feature extents.  For example "FDb.InfiniumMethylation.hg19:get27k"
is valid.  Note that in this case, gchr must have prefix "chr".
}
  \item{snpannopk}{
character string naming a Bioconductor SNPlocs.* SNP metadata package
}
  \item{as.GRangesList}{
logical telling whether a GRangesList should be returned
}
  \item{excludeRadius}{
numeric or NULL: radius of interval around gene extent from
which SNP will be excluded, required to be less than \code{radius}
}
}
\details{
This is a utility that the developer would rather
not export.  The complexity of harmonizing queries among
probe and SNP annotation resources leads to a somewhat fragile
product.  Users who have their own gene ranges and SNP
locations can examine the namelist component of the
output of the
default call to see what is expected for
the *.cis.eQTLs function.  For the set of
chromosomes to be analyzed, there will be a list
of chromosome specific namelist-like lists. 
}
\value{
Instance of \code{cisMap} class, which will retain
SNP location, gene range, and radius information for
auditing.
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
 \dontrun{
  getCisMap()
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
