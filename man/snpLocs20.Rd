\name{snpLocs20}
\alias{snpLocs20}
\docType{data}
\title{
prototype SNP location instance for use with GGtools
}
\description{
prototype SNP location instance for use with GGtools
}
\usage{data(snpLocs20)}
\format{
  The format is:
Formal class 'UCSCData' [package "rtracklayer"] with 6 slots
  ..@ trackLine      :Formal class 'BasicTrackLine' [package "rtracklayer"] with 12 slots \cr
  .. .. ..@ itemRgb    : logi(0)  \cr
  .. .. ..@ useScore   : logi(0)  \cr
  .. .. ..@ group      : chr(0)  \cr
  .. .. ..@ db         : chr(0)  \cr
  .. .. ..@ offset     : num(0)  \cr
  .. .. ..@ url        : Named chr " " \cr
  .. .. .. ..- attr(*, "names")= chr "url" \cr
  .. .. ..@ htmlUrl    : chr(0)  \cr
  .. .. ..@ name       : Named chr "snp130" \cr
  .. .. .. ..- attr(*, "names")= chr "name" \cr
  .. .. ..@ description: Named chr "snp130" \cr
  .. .. .. ..- attr(*, "names")= chr "description" \cr
  .. .. ..@ visibility : Named chr "1" \cr
  .. .. .. ..- attr(*, "names")= chr "visibility" \cr
  .. .. ..@ color      : int(0)  \cr
  .. .. ..@ priority   : num(0)  \cr
  ..@ ranges         :Formal class 'CompressedIRangesList' [package "IRanges"] with 5 slots \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ elementType    : chr "IRanges" \cr
  .. .. ..@ metadata       :List of 1 \cr
  .. .. .. ..$ universe: chr "hg18" \cr
  .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots \cr
  .. .. .. .. ..@ end            : int 450693 \cr
  .. .. .. .. ..@ NAMES          : chr "chr20" \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ unlistData     :Formal class 'IRanges' [package "IRanges"] with 6 slots \cr
  .. .. .. .. ..@ start          : int [1:450693] 60492 60572 60646 60705 61098 61605 61795 62100 62291 62731 ... \cr
  .. .. .. .. ..@ width          : int [1:450693] 1 1 1 0 1 1 1 1 0 1 ... \cr
  .. .. .. .. ..@ NAMES          : NULL \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. ..@ metadata       : list() \cr
  ..@ values         :Formal class 'CompressedSplitDataFrameList' [package "IRanges"] with 5 slots \cr
  .. .. ..@ elementMetadata: NULL \cr
  .. .. ..@ elementType    : chr "DataFrame" \cr
  .. .. ..@ metadata       : list() \cr
  .. .. ..@ partitioning   :Formal class 'PartitioningByEnd' [package "IRanges"] with 5 slots \cr
  .. .. .. .. ..@ end            : int 450693 \cr
  .. .. .. .. ..@ NAMES          : chr "chr20" \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ elementType    : chr "integer" \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. ..@ unlistData     :Formal class 'DataFrame' [package "IRanges"] with 6 slots \cr
  .. .. .. .. ..@ rownames       : NULL \cr
  .. .. .. .. ..@ nrows          : int 450693 \cr
  .. .. .. .. ..@ elementMetadata: NULL \cr
  .. .. .. .. ..@ elementType    : chr "ANY" \cr
  .. .. .. .. ..@ metadata       : list() \cr
  .. .. .. .. ..@ listData       :List of 3 \cr
  .. .. .. .. .. ..$ name  : chr [1:450693] "rs35078228" "rs28753379" "rs28579812" "rs35616340" ... \cr
  .. .. .. .. .. ..$ score : num [1:450693] 0 0 0 0 0 0 0 0 0 0 ... \cr
  .. .. .. .. .. ..$ strand: chr [1:450693] "+" "+" "+" "+" ... \cr
  ..@ elementMetadata: NULL \cr
  ..@ elementType    : chr "ANY" \cr
  ..@ metadata       : list() \cr
}
\details{
derived from UCSC table for snp130
}
\source{
snp130 table in hg19 UCSC table set
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(snpLocs20)
snpLocs20
}
\keyword{datasets}
