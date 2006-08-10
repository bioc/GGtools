\name{HM2rac}
\alias{HM2rac}
\title{compute rare allele count from a hapmap file}
\description{compute rare allele count from a hapmap file}
\usage{
HM2rac(fn, comment.char="", kppref="NA")
}
\arguments{
  \item{fn}{ filename string identifying hapmap genotyping file}
  \item{comment.char}{ string identifying comment character
token for read.table, should not be sharp sign because hapmap
will use this in the header line }
  \item{kppref}{ a character token that must be present to
identify a sample header column.  I.e., if kppref is "NA"
then a column with header tag NAnnnnn will be kept but
a column with header chrom will be dropped.
}
\value{named character matrix with one row per snp and one column
per sample, with snps named by rs numbers.
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
fi = system.file("fileDemos/geno20demo.txt")
HM2rac(fi)[1:5,1:7]
}
\keyword{ models }
