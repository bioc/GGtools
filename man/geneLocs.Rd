\name{geneLocs}
%c15GGceu.rda    c7GGceu.rda     chr7meta.rda    racHMceu20.rda
%c20GG.rda       chr15meta.rda   geneLocs.rda    racHMceu7.rda
%c20GGceu.rda    chr20meta.rda   racHMceu15.rda
%
\alias{geneLocs}
\alias{geneLocs_hsa}
\alias{geneLocs_mmu}
\title{gene metadata from NCBI}
\description{
gene metadata from NCBI
}
\usage{
data(geneLocs)
data(geneLocs_hsa)
data(geneLocs_mmu)
}
\details{
A \code{\link{data.frame}} instance with information relating
gene symbol to chromosome and providing start and end positional
information.
}
\value{
instance of class \code{\link{data.frame}}
}
\references{ 
NCBI
\url{ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/mapview/seq_gene.md.gz}
for the human data.

\url{ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/mapview/seq_gene.md.gz}
for the mouse data.
 }
\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{The geneLocs structure should be available
for multiple species.  Thus the three letter code
suffix will be present on the standard offering in the future.  

There is absolutely no guarantee that the information in these
structures is at all accurate.  The provenance subdirectory of the
package includes notes on how the files were derived.
The package author finds it difficult to obtain clear information on the
state of the art on gene location assertions; suggestions of more
accurate sources are welcome.
}


%\seealso{  }

\examples{
data(geneLocs_hsa)
# be patient
geneLocs[1:4,]
}
\keyword{ models }
