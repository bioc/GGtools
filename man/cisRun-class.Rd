\name{cisRun-class}
\Rdversion{1.1}
\docType{class}
\alias{cisRun-class}

\title{Class \code{"cisRun"}}
\description{
manage results of All.cis eQTL analysis
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cisRun", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{seqnames}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{ranges}:}{Object of class \code{"IRanges"} ~~ }
    \item{\code{strand}:}{Object of class \code{"Rle"} ~~ }
    \item{\code{elementMetadata}:}{Object of class \code{"DataFrame"} ~~ }
    \item{\code{seqinfo}:}{Object of class \code{"Seqinfo"} ~~ }
    \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{GRanges}"}, directly.
Class \code{"\linkS4class{GenomicRanges}"}, by class "GRanges", distance 2.
Class \code{"\linkS4class{Vector}"}, by class "GRanges", distance 3.
Class \code{"\linkS4class{GenomicRangesORmissing}"}, by class "GRanges", distance 3.
Class \code{"\linkS4class{GenomicRangesORGRangesList}"}, by class "GRanges", distance 3.
%Class \code{"\linkS4class{RangedDataORGenomicRanges}"}, by class "GRanges", distance 3.
Class \code{"\linkS4class{Annotated}"}, by class "GRanges", distance 4.
}
\section{Methods}{
No methods defined with class "cisRun" in the signature.
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
intent is to simplify output of cis eQTL testing in a GRanges
instance
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("cisRun")
}
\keyword{classes}