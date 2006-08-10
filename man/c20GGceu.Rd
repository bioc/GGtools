\name{c20GGceu}
%c15GGceu.rda    c7GGceu.rda     chr7meta.rda    racHMceu20.rda
%c20GG.rda       chr15meta.rda   geneLocs.rda    racHMceu7.rda
%c20GGceu.rda    chr20meta.rda   racHMceu15.rda
%
\alias{c20GGceu}
\alias{c7GGceu}
\alias{c15GGceu}
\title{representations of HapMap snp data + expression data}
\description{
representations of HapMap snp data + expression data
}
\usage{
data(c20GGceu)
}
\details{
These instances of class \code{\link{ggExprSet}} are
created from two basic sources.  

First, V. Cheung and R. Spielman
provided images of the expression data for the 58 CEPH individuals
studied in their Oct 2005 Nature letter.  This is expression
data created using the hgfocus platform.  Details of preprocessing
or provenance still to be supplied; this particular data structure
was distributed at CDATA-06 in June 2006 at Cold Spring Harbod.

Second, data on forward non-redundant SNPs in 48 individuals
among the 58 mentioned above were obtained from the HapMap
ftp site in late July 2006.  Full provenance information
still to be supplied.
}
\value{
instances of class \code{\link{ggExprSet}}
}
\references{ 
Cheung VG., Spielman RS., Ewens KG., Weber TM., Morley M & Burdick JT.: Mapping determinants of human gene expression by regional and whole genome association. Nature, 437: 1365-1369, 2005
 }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(c20GGceu)
# be patient
exprs(c20GGceu)[1:4,1:4]
snps(c20GGceu)[1:4,1:4]
}
\keyword{ models }
