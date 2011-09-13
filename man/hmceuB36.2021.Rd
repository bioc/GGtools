\name{hmceuB36.2021}
\alias{ex}
\docType{data}
\title{ two chromosomes of genotype data and full expression data for CEPH CEU hapmap data}
\description{ two chromosomes of genotype data and full expression data for CEPH CEU hapmap data}
\usage{
data(ex)  # not intended for direct use, instead use getSS
# getSS("GGtools", c("20", "21"))
}
\format{
  The format is:
Formal class 'smlSet' [package "GGBase"] with 9 slots

  ..@ smlEnv           :<environment: 0x3902e98> 

  ..@ annotation       : chr "illuminaHumanv1.db"

  ..@ chromInds        : num [1:2] 20 21

  ..@ organism         : chr "Hs"

  ..@ assayData        :<environment: 0x3c96504> 

  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots

  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots

  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots

  ..@ ...classVersion..:Formal class 'Versions' [package "Biobase"] with 1 slots
}
\details{
Until Sept 2011 this object was serialized as an smlSet instance.
Now the ExpressionSet component is serialized in data (as eset.rda, with
contents ex), and the genotype data are in inst/parts as SnpMatrix instances.
}
%\source{
%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
if (!exists("hmceuB36.2021")) hmceuB36.2021 <- getSS("GGtools", c("20", "21"))
validObject(hmceuB36.2021)
}
\keyword{datasets}
