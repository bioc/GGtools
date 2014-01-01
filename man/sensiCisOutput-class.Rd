\name{sensiCisOutput-class}
\Rdversion{1.1}
\docType{class}
\alias{sensiCisOutput-class}
\alias{show,sensiCisOutput-method}

\title{Class \code{"sensiCisOutput"}}
\description{
This class helps to manage the results from a collection
of cis-eQTL searches.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sensiCisOutput", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{byGene}:}{Object of class \code{"GRanges"}, organized
to provide ranges for genes and their best associated
cis SNP  }
    \item{\code{bySNP}:}{Object of class \code{"GRanges"} organized
to provide easy access to genomic coordinates of SNP found to
be most strongly associated with a gene in cis }
    \item{\code{tabAtFDRB}:}{Object of class \code{"ANY"} a flattened
table that defines tuning parameters and eQTL yield for a collection
of searches }
    \item{\code{input}:}{Object of class \code{"sensiCisInput"} :
object that describes the files and parameter settings used
for the sensitivity analysis }
    \item{\code{thecall}:}{Object of class \code{"call"}: the
call generating this instance }
    \item{\code{fdrbound}:}{Object of class \code{"numeric"}: gives
the upper bound on FDR for declaring an eQTL  }
    \item{\code{sessionInfo}:}{Object of class \code{"ANY"}: describes
state of system in which the object was made.  }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "sensiCisOutput")}: concise
rendering with hints }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("sensiCisOutput")
}
\keyword{classes}
