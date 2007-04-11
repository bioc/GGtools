\name{buildOGTES}
\alias{buildOGTES}
\title{ combine Expression and SNP data into oGtypeExSet instance }
\description{
 combine Expression and SNP data into oGtypeExSet instance 
}
\usage{
buildOGTES(exSet, ...)
}
\arguments{
  \item{exSet}{ instance of ExpressionSet }
  \item{\dots}{ named sequence of SnpCallSet instances; suggest names Nsp, Sty, etc. }
}
\details{
combines phenoData using combine method
}
\value{
instance of oGtypeExSet
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
getClass("oGtypeExSet")
# if you have Jeff's nspCalls, styCalls, x3pSet, then
# buildOGTES( x3pSet, Nsp=nspCalls, Sty=styCalls )
# will give a nice show
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
