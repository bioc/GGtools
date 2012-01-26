\name{ex}
\alias{ex}
\docType{data}
\title{
ExpressionSet instance for illustrating integrative smlSet container
}
\description{
ExpressionSet instance for illustrating integrative smlSet container
}
\usage{data(ex)}
\format{
  The format is:
Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots \cr
  .. .. ..@ name             : chr "" \cr
  .. .. ..@ lab              : chr "" \cr
  .. .. ..@ contact          : chr "" \cr
  .. .. ..@ title            : chr "" \cr
  .. .. ..@ abstract         : chr "" \cr
  .. .. ..@ url              : chr "" \cr
  .. .. ..@ pubMedIds        : chr "" \cr
  .. .. ..@ samples          : list() \cr
  .. .. ..@ hybridizations   : list() \cr
  .. .. ..@ normControls     : list() \cr
  .. .. ..@ preprocessing    : list() \cr
  .. .. ..@ other            : list() \cr
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots \cr
  .. .. .. .. ..@ .Data:List of 2 \cr
  .. .. .. .. .. ..$ : int [1:3] 1 0 0 \cr
  .. .. .. .. .. ..$ : int [1:3] 1 1 0 \cr
  ..@ assayData        :<environment: 0x10bf12948>  \cr
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots \cr
  .. .. ..@ varMetadata      :'data.frame':	7 obs. of  1 variable: \cr
  .. .. .. ..$ labelDescription: chr [1:7] "hapmap family id" "hapmap person id" "id of mother of this person" "id of father of this person" ... \cr
  .. .. ..@ data             :'data.frame':	90 obs. of  7 variables: \cr
  .. .. .. ..$ famid    : int [1:90] 1341 1341 1341 1340 1340 1340 1340 1340 1341 1341 ... \cr
  .. .. .. ..$ persid   : int [1:90] 14 2 13 9 10 2 11 1 11 1 ... \cr
  .. .. .. ..$ mothid   : int [1:90] 0 14 0 0 0 12 0 10 0 12 ... \cr
  .. .. .. ..$ fathid   : int [1:90] 0 13 0 0 0 11 0 9 0 11 ... \cr
  .. .. .. ..$ sampid   : Factor w/ 90 levels "NA06985","NA06991",..: 1 2 3 4 5 6 7 8 9 10 ... \cr
  .. .. .. ..$ isFounder: logi [1:90] TRUE FALSE TRUE TRUE TRUE FALSE ... \cr
  .. .. .. ..$ male     : logi [1:90] FALSE FALSE TRUE TRUE FALSE FALSE ... \cr
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns" \cr
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots \cr
  .. .. .. .. ..@ .Data:List of 1 \cr
  .. .. .. .. .. ..$ : int [1:3] 1 1 0 \cr
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots \cr
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable: \cr
  .. .. .. ..$ labelDescription: chr(0)  \cr
  .. .. ..@ data             :'data.frame':	47293 obs. of  0 variables \cr
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns" \cr
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots \cr
  .. .. .. .. ..@ .Data:List of 1 \cr
  .. .. .. .. .. ..$ : int [1:3] 1 1 0 \cr
  ..@ annotation       : chr "illuminaHumanv1.db" \cr
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots \cr
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable: \cr
  .. .. .. ..$ labelDescription: chr(0)  \cr
  .. .. ..@ data             :'data.frame':	90 obs. of  0 variables \cr
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns" \cr
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots \cr
  .. .. .. .. ..@ .Data:List of 1 \cr
  .. .. .. .. .. ..$ : int [1:3] 1 1 0 \cr
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots \cr
  .. .. ..@ .Data:List of 4 \cr
  .. .. .. ..$ : int [1:3] 2 14 0 \cr
  .. .. .. ..$ : int [1:3] 2 13 7 \cr
  .. .. .. ..$ : int [1:3] 1 3 0 \cr
  .. .. .. ..$ : int [1:3] 1 0 0 \cr
}
\details{
Expression data harvested in 2007 from GENEVAR

ftp://ftp.sanger.ac.uk/pub/genevar/CEU_parents_norm_march2007.zip

}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(ex)
}
\keyword{datasets}
