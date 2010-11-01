\name{readGVF}
\alias{readGVF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read and model a GVF (Genome Variation Format) file (see Reese M et al, Genome Biology 2010 11:R88)
}
\description{
read and model a GVF (Genome Variation Format) file (see Reese M et al, Genome Biology 2010 11:R88)
}
\usage{
readGVF(gzfname, ncheck = 100, ndata = NULL, els = c("ID", "Variant_seq", "Reference_seq", "Total_reads", "Genotype"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gzfname}{
path to a gzipped text file
}
  \item{ncheck}{
guess at the number of records consumed by pragmas and comments at the head of the file
}
  \item{ndata}{
number of data lines to be processed; if NULL, all
}
  \item{els}{
key names for the 9th field
}
}
\details{
small footprint, slow parsing into GRanges instance
}
\value{
\code{\link[GenomicRanges]{GRanges-class}} instance
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Vince Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
x = readGVF( system.file("GVF/10Gen_NA19240_SNVexc.gvf.gz", package="GGtools") )
x
}
\keyword{ models }
