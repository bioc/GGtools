\name{Strains2rac}
\alias{Strains2rac}
\alias{INBREDSworkflow}
\alias{INB34snpMeta}
\title{convert a Wellcome 'Strains' genotyping file to rare allele count form}
\description{convert a Wellcome 'Strains' genotyping file to rare allele count form,
or create a racExSet using such a file
}
\usage{
Strains2rac(strfile)
INBREDSworkflow(inbfile, emat, estrains, pd, mi, anno, fixup= NULL, fixchr =
     function(x) gsub("_random", "", x))
}
\arguments{
  \item{strfile}{ name of the text file to be converted}
  \item{inbfile}{ name of the text file to be converted}
  \item{emat}{ matrix of expression data }
  \item{estrains}{ character vector of strains corresponding to
columns of \code{emat}}
  \item{pd}{ \code{\link[Biobase]{AnnotatedDataFrame-class}} instance for \code{emat}}
  \item{mi}{ \code{\link[Biobase]{MIAME-class}} instance}
  \item{anno}{ character string for annotation environment package}
  \item{fixup}{ transformation from columns of strain genotype file \code{strfile}
 to \code{estrains} tokens, if necessary}
  \item{fixchr}{ modify the chromosome tags if necessary }
}
\details{
}
\value{
\code{Strains2rac} creates a list with genotype matrix, chr and pos

\code{INBREDSworkflow} creates a list with components \code{\link{racExSet}} and
\code{\link{snpMetaWhole}} (because the INBREDS files are whole-genome, not
chromosome specific)
}
\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{See the commented example for build 34 Wellcome Inbreds file for
BXD strains, \code{INB34snpMeta}.  }


\seealso{  }

\examples{
fi = system.file("fileDemos/StrainInit.txt", package="GGtools")
dem = Strains2rac(fi)
names(dem)
dim(dem[[1]])
data(gse2031GG)
fixbxd = function(x) {
# because such various labels are used in the INBREDS file
 gsub("BXD\\\.", "BXD", gsub(".Ty", "", gsub(".TyJ", "", x)))
}
INBREDSworkflow( fi, exprs(gse2031GG), as.character(gse2031GG$str),
	phenoData(gse2031GG), new("MIAME"), "mgu74av2", fixbxd )
#
# the real INB34snpMeta structure is built by
# fixchr = function(x) gsub("_random", "", x)
#INBREDSworkflow( "Strains-0.09052005.txt", exprs(gse2031GG), as.character(gse2031GG$str),
#	phenoData(gse2031GG), new("MIAME"), "mgu74av2", fixbxd, fixchr )
}
\keyword{ models }
