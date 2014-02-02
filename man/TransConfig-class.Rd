\name{TransConfig-class}
\Rdversion{1.1}
\docType{class}
\alias{TransConfig-class}
\alias{batchsize,TransConfig-method}
\alias{batchsize<-,TransConfig,integer-method}
\alias{gbufsize,TransConfig-method}
\alias{gbufsize<-,TransConfig,integer-method}
\alias{show,TransConfig-method}
\alias{snpchr,TransConfig-method}
\alias{snpchr<-,TransConfig,character-method}

\title{Class \code{"TransConfig"}}
\description{
Instances from this class can be input to the transScores
function to control a trans-eQTL search.
}
\section{Objects from the Class}{
Instances from this class can be input to the transScores
function to control a trans-eQTL search.
Objects can be created by calls of the form \code{new("TransConfig")}.
}
\section{Slots}{
  \describe{
    \item{\code{snpchr}:}{Object of class \code{"character"} 
identifies the name of the chromosome harboring
SNP that will all be used (subject
to filtering by smFilter function) in transcriptome-wide searches
for associated transcripts }
    \item{\code{gbufsize}:}{Object of class \code{"integer"},
scores for the top gbufsize genes are retained as the search proceeds  }
    \item{\code{batchsize}:}{Object of class \code{"integer"} 
used in processing ff-based archives for association scores }
    \item{\code{smpack}:}{Object of class \code{"character"},
tells the name of the installed package used for retrieval
of expression-genotype data using \code{\link[GGBase]{getSS}}  }
    \item{\code{rhs}:}{Object of class \code{"formula"}, formula
used in \code{\link[snpStats]{snp.rhs.tests}}; typically not used.
If plug-in FDR is desired,
adjustments should be executed in a regressOut call .}
%    \item{\code{nperm}:}{Object of class \code{"integer"} ~~ }
    \item{\code{folderStem}:}{Object of class \code{"character"},
name of a folder where interim results are sequestered.  }
    \item{\code{radius}:}{Object of class \code{"integer"},
defines region around SNP within which genes are considered 'cis' so
tests are not conducted. }
    \item{\code{shortfac}:}{Object of class \code{"integer"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{chrnames}:}{Object of class \code{"character"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{smchrpref}:}{Object of class \code{"character"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{gchrpref}:}{Object of class \code{"character"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{schrpref}:}{Object of class \code{"character"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{geneApply}:}{Object of class \code{"function"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{geneannopk}:}{Object of class \code{"character"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{snpannopk}:}{Object of class \code{"character"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{smFilter}:}{Object of class \code{"function"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{exFilter}:}{Object of class \code{"function"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{keepMapCache}:}{Object of class \code{"logical"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{SSgen}:}{Object of class \code{"function"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{excludeRadius}:}{Object of class \code{"integerOrNULL"} see
documentation for \code{\link{CisConfig-class}} }
    \item{\code{estimates}:}{Object of class \code{"logical"} see
documentation for \code{\link{CisConfig-class}} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{CisConfig}"}, directly.
}
\section{Methods}{
  \describe{
    \item{batchsize}{\code{signature(x = "TransConfig")}: ... }
    \item{batchsize<-}{\code{signature(object = "TransConfig", value = "integer")}: ... }
    \item{gbufsize}{\code{signature(x = "TransConfig")}: ... }
    \item{gbufsize<-}{\code{signature(object = "TransConfig", value = "integer")}: ... }
    \item{show}{\code{signature(object = "TransConfig")}: ... }
    \item{snpchr}{\code{signature(x = "TransConfig")}: ... }
    \item{snpchr<-}{\code{signature(object = "TransConfig", value = "character")}: ... }
	 }
}
%\references{
%%%  ~~put references to the literature/web site here~~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("TransConfig")
}
\keyword{classes}
