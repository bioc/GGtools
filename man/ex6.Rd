\name{ex6}
\alias{ex6}
\docType{data}
\title{
example exon region data
}
\description{
example exon region data
}
\usage{data(ex6)}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 7 slots
  ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 5 slots
  .. .. ..@ values         : Factor w/ 49 levels "chr1","chr1_random",..: 36
  .. .. ..@ lengths        : int 12974
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:12974] 237101 249628 256880 280114 290854 293103 293769 293769 295822 336752 ...
  .. .. ..@ width          : int [1:12974] 460 34 83 50 75 172 73 2585 534 58 ...
  .. .. ..@ NAMES          : NULL
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ elementType    : chr "integer"
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "IRanges"] with 5 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 1 2
  .. .. ..@ lengths        : int [1:2] 6235 6739
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ metadata       : list()
  ..@ seqlengths     : Named int [1:49] 247249719 1663265 135374737 113275 134452384 215294 132349534 114142980 186858 106368585 ...
  .. ..- attr(*, "names")= chr [1:49] "chr1" "chr1_random" "chr10" "chr10_random" ...
  ..@ elementMetadata:Formal class 'DataFrame' [package "IRanges"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 12974
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ metadata       : list()
  .. .. ..@ listData       :List of 1
  .. .. .. ..$ exon_id: int [1:12974] 81518 81519 81520 81521 81522 81523 81524 81526 81525 81527 ...
  ..@ elementType    : chr "ANY"
  ..@ metadata       : list()
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ex6)
## maybe str(ex6) ; plot(ex6) ...
}
\keyword{datasets}
