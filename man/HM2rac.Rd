\name{HM2rac}
\alias{HM2rac}
\alias{thinHM2rac}
\alias{thinHM2meta}
\alias{df2snpMeta}
\alias{getRare}
\alias{countRare}
\title{compute rare allele count from a hapmap file}
\description{compute rare allele count from a hapmap file}
\usage{
HM2rac(fn, comment.char="", kppref="NA")
thinHM2rac(gzfn)
countRare(x)
getRare(x)
}
\arguments{
  \item{fn}{ filename string identifying hapmap genotyping file}
  \item{gzfn}{ filename string identifying gzipped hapmap genotyping file}
  \item{x}{ string vector of simple genotype calls (e.g., "AA") }
  \item{comment.char}{ string identifying comment character
token for read.table, should not be sharp sign because hapmap
will use this in the header line }
  \item{kppref}{ a character token that must be present to
identify a sample header column.  I.e., if kppref is "NA"
then a column with header tag NAnnnnn will be kept but
a column with header chrom will be dropped.}
}
\value{HM2rac returns a list including elements
raremat: named character matrix with one row per snp and one column
per sample, with snps named by rs numbers

alleles: named vector showing existing alleles at snp

rarealleles: named vector identifying base of rare allele

anno: character string for annotation indirection

countRare and getRare are simple string processing functions

\code{thinHM2meta} will create a \code{\link{snpMeta-class}} structure from a HapMap gzipped genotyping
file.
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
fi = system.file("fileDemos/geno20demo.txt.gz", package="GGtools")
racstr =  thinHM2rac(fi)
names(racstr)
racstr$rarecount[1:5,1:7]
racstr$rarebase[1:6]
racstr$alleles[1:6]
getRare(c("AA", "AG", "GG", "GG"))
countRare(c("AA", "AG", "GG", "GG"))
}
\keyword{ models }
