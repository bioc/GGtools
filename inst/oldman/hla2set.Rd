\name{hla2set}
\alias{hla2set}
\docType{data}
\title{ a gene set of 9 genes from human HLA2 locus }
\description{ a gene set of 9 genes from human HLA2 locus }

\usage{data(hla2set)}
\format{
  The format is:
Formal class 'GeneSet' [package "GSEABase"] with 13 slots

  ..@ geneIdType      :Formal class 'SymbolIdentifier' [package "GSEABase"] with 2 slots

  .. .. ..@ type      :Formal class 'ScalarCharacter' [package "Biobase"] with 1 slots

and so on.

See \code{\link[GSEABase]{GeneSet-class}} for additional information.
}
\details{
This set of 9 genes related to human HLA2 locus was used in the
2009 Bioinformatics Application Note by Carey, Davis et al.
}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
if ("GSEABase" \%in\% installed.packages()[,1]) {
 load(system.file("genesets/hla2set.rda", package="GGtools"))
 hla2set
 }
}
\keyword{datasets}
