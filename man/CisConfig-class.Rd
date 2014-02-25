\name{CisConfig-class}
\Rdversion{1.1}
\docType{class}
\alias{CisConfig-class}
\alias{genome,CisConfig-method}
\alias{genome<-,CisConfig,character-method}
\alias{chrnames,CisConfig-method}
\alias{chrnames<-,CisConfig,character-method}
\alias{estimates,CisConfig-method}
\alias{estimates<-,CisConfig,logical-method}
\alias{excludeRadius,CisConfig-method}
\alias{excludeRadius<-,CisConfig,integer-method}
\alias{exFilter,CisConfig-method}
\alias{exFilter<-,CisConfig,function-method}
\alias{gchrpref,CisConfig-method}
\alias{gchrpref<-,CisConfig,character-method}
\alias{geneannopk,CisConfig-method}
\alias{geneannopk<-,CisConfig,character-method}
\alias{geneApply,CisConfig-method}
\alias{geneApply<-,CisConfig,function-method}
\alias{initialize,CisConfig-method}
\alias{keepMapCache,CisConfig-method}
\alias{keepMapCache<-,CisConfig,logical-method}
\alias{radius,CisConfig-method}
\alias{radius<-,CisConfig,integer-method}
\alias{rhs,CisConfig-method}
\alias{rhs<-,CisConfig,function-method}
\alias{schrpref,CisConfig-method}
\alias{schrpref<-,CisConfig,character-method}
\alias{shortfac,CisConfig-method}
\alias{shortfac<-,CisConfig,integer-method}
\alias{show,CisConfig-method}
\alias{smchrpref,CisConfig-method}
\alias{smchrpref<-,CisConfig,character-method}
\alias{smFilter,CisConfig-method}
\alias{smFilter<-,CisConfig,function-method}
\alias{snpannopk,CisConfig-method}
\alias{snpannopk<-,CisConfig,character-method}
\alias{SSgen,CisConfig-method}
\alias{SSgen<-,CisConfig,function-method}
\alias{chrnames}
\alias{chrnames<-}
\alias{estimates}
\alias{estimates<-}
\alias{excludeRadius}
\alias{excludeRadius<-}
\alias{exFilter}
\alias{exFilter<-}
\alias{gchrpref}
\alias{gchrpref<-}
\alias{geneannopk}
\alias{geneannopk<-}
\alias{geneApply}
\alias{geneApply<-}
\alias{initialize}
\alias{keepMapCache}
\alias{keepMapCache<-}
\alias{radius}
\alias{radius<-}
\alias{rhs}
\alias{rhs<-}
\alias{schrpref}
\alias{schrpref<-}
\alias{shortfac}
\alias{shortfac<-}
\alias{show}
\alias{smchrpref}
\alias{smchrpref<-}
\alias{smFilter}
\alias{smFilter<-}
\alias{snpannopk}
\alias{snpannopk<-}
\alias{SSgen}
\alias{SSgen<-}
\alias{genome}
\alias{genome<-}
\alias{extraProps}
\alias{extraProps<-}
\alias{extraProps,CisConfig-method}
\alias{extraProps<-,CisConfig,function-method}

\alias{excludeRadius<-,CisConfig,integerOrNULL-method}
\alias{folderStem,CisConfig-method}
\alias{folderStem<-,CisConfig,character-method}
\alias{nperm,CisConfig-method}
\alias{nperm<-,CisConfig,integer-method}
\alias{rhs<-,CisConfig,formula-method}
\alias{smpack,CisConfig-method}
\alias{smpack<-,CisConfig,character-method}

\alias{excludeRadius<-}
\alias{folderStem}
\alias{folderStem<-}
\alias{nperm}
\alias{nperm<-}
\alias{rhs<-}
\alias{smpack}
\alias{smpack<-}


\title{Class \code{"CisConfig"}}
\description{
Object specifying configuration of cis-eQTL search, to be used
with All.cis
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CisConfig")}.
Use replacement methods to update the fields.
}
\section{Slots}{
  \describe{
    \item{\code{smpack}:}{character string identifying package
holding the expression and genotype data; see \code{\link[GGBase]{getSS}}}
    \item{\code{genome}:}{character string identifying genome build in use}
    \item{\code{rhs}:}{Object of class \code{"formula"} right hand
side for calls to \code{\link[snpStats]{snp.rhs.tests}} }
    \item{\code{nperm}:}{Object of class \code{"integer"} 
number of permutations for plug in FDR}
    \item{\code{folderStem}:}{Object of class \code{"character"} 
string used for scratch space folders, relative to current folder }
    \item{\code{radius}:}{Object of class \code{"integer"} radius
of search }
    \item{\code{shortfac}:}{Object of class \code{"integer"} 
scores are scaled up by this factor so that precision can
be retained in short integer representation }
    \item{\code{chrnames}:}{Object of class \code{"character"} 
string identifying chromosome label used in gene annotation
retrieval -- typically length 1 }
    \item{\code{smchrpref}:}{Object of class \code{"character"}
prefix to be attached to chromosome label in \code{chrnames}
to pick out the element of \code{\link[GGBase]{smlSet-class}} instance
used in testing }
    \item{\code{gchrpref}:}{Object of class \code{"character"} prefix
on chrnames token to be used for gene location retrievals }
    \item{\code{schrpref}:}{Object of class \code{"character"} prefix
on chrnames token to be used with SNPlocs package for
retrieval of SNP locations }
    \item{\code{geneApply}:}{Object of class \code{"function"} iterator
over genes, could be lapply or mclapply }
    \item{\code{geneannopk}:}{Object of class \code{"character"} Bioconductor
annotation package for gene locations, typically for expression array }
    \item{\code{snpannopk}:}{Object of class \code{"character"} Bioconductor
dbSNP annotation package }
    \item{\code{smFilter}:}{Object of class \code{"function"} function
to be applied to smlSet instance that yields an smlSet instance with
required contents; could apply MAF restriction for example by calling
MAFfilter }
    \item{\code{exFilter}:}{Object of class \code{"function"} function
that is run right after smlSet is materialized, permitting replacement
or filtering of expression data, when, for example, the ExpressionSet
includes multiple tissue types }
    \item{\code{keepMapCache}:}{Object of class \code{"logical"} for
enhancing processing of gene-SNP cis mapping with a global cache }
    \item{\code{SSgen}:}{Object of class \code{"function"} function
that accepts name of an expression+SnpMatrix package (as generated
by \code{\link{externalize}}), a chromosome tag (\code{chrnames}
prefixed by \code{smchrpref}), and a function, and returns
an smlSet instance }
    \item{\code{excludeRadius}:}{Object of class \code{"integerOrNULL"}
which will determine what interval about the gene is excluded for
cis testing; 0 should exclude all within-gene SNP, but needs
testing  }
    \item{\code{estimates}:}{Object of class \code{"logical"} if TRUE,
estimates and standard errors (expanded and reduced in storage
as a short int, using \code{shortfac}) are generated and retained }
    \item{\code{extraProps}:}{Object of class \code{"function"} 
this function is applied to the cisScores output
before it is returned, to bind additional metadata
to the ranges if desired.  Defaults to function(x)x. }
    \item{\code{useME}:}{Object of class \code{"logical"} if TRUE,
use the statistics generated by \code{\link[MatrixEQTL]{Matrix_eQTL_engine}} for
association testing.}
    \item{\code{MEpvot}:}{Object of class \code{"numeric"} used
if \code{useME} slot is set to TRUE: p-value output threshold
for retaining association test statistic generated by
\code{\link[MatrixEQTL]{Matrix_eQTL_engine}}; defaults to 0.5.  Higher values
lead to higher volumes and longer times to completion.}
  }
}
\section{Methods}{
  \describe{
    \item{chrnames}{\code{signature(x = "CisConfig")}: ... }
    \item{chrnames<-}{\code{signature(object = "CisConfig", value = "character")}: ... }
    \item{estimates}{\code{signature(x = "CisConfig")}: ... }
    \item{estimates<-}{\code{signature(object = "CisConfig", value = "logical")}: ... }
    \item{excludeRadius}{\code{signature(x = "CisConfig")}: ... }
    \item{excludeRadius<-}{\code{signature(object = "CisConfig", value = "integer")}: ... }
    \item{exFilter}{\code{signature(x = "CisConfig")}: ... }
    \item{exFilter<-}{\code{signature(object = "CisConfig", value = "function")}: ... }
    \item{gchrpref}{\code{signature(x = "CisConfig")}: ... }
    \item{gchrpref<-}{\code{signature(object = "CisConfig", value = "character")}: ... }
    \item{geneannopk}{\code{signature(x = "CisConfig")}: ... }
    \item{geneannopk<-}{\code{signature(object = "CisConfig", value = "character")}: ... }
    \item{geneApply}{\code{signature(x = "CisConfig")}: ... }
    \item{geneApply<-}{\code{signature(object = "CisConfig", value = "function")}: ... }
    \item{initialize}{\code{signature(.Object = "CisConfig")}: ... }
    \item{keepMapCache}{\code{signature(x = "CisConfig")}: ... }
    \item{keepMapCache<-}{\code{signature(object = "CisConfig", value = "logical")}: ... }
    \item{radius}{\code{signature(x = "CisConfig")}: ... }
    \item{radius<-}{\code{signature(object = "CisConfig", value = "integer")}: ... }
    \item{rhs}{\code{signature(x = "CisConfig")}: ... }
    \item{rhs<-}{\code{signature(object = "CisConfig", value = "function")}: ... }
    \item{schrpref}{\code{signature(x = "CisConfig")}: ... }
    \item{schrpref<-}{\code{signature(object = "CisConfig", value = "character")}: ... }
    \item{shortfac}{\code{signature(x = "CisConfig")}: ... }
    \item{shortfac<-}{\code{signature(object = "CisConfig", value = "integer")}: ... }
    \item{show}{\code{signature(object = "CisConfig")}: ... }
    \item{smchrpref}{\code{signature(x = "CisConfig")}: ... }
    \item{smchrpref<-}{\code{signature(object = "CisConfig", value = "character")}: ... }
    \item{smFilter}{\code{signature(x = "CisConfig")}: ... }
    \item{smFilter<-}{\code{signature(object = "CisConfig", value = "function")}: ... }
    \item{snpannopk}{\code{signature(x = "CisConfig")}: ... }
    \item{snpannopk<-}{\code{signature(object = "CisConfig", value = "character")}: ... }
    \item{SSgen}{\code{signature(x = "CisConfig")}: ... }
    \item{SSgen<-}{\code{signature(object = "CisConfig", value = "function")}: ... }
	 }
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("CisConfig")
}
\keyword{classes}
