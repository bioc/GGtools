\name{supp2SQLite}
\alias{supp2SQLite}
\title{function to transform snp.support data
from snpMatrix read.HapMap.data into SQLite table  }
\description{function to transform snp.support data
from snpMatrix read.HapMap.data into SQLite table  }
\usage{
supp2SQLite(df, tabname, fn, ...)
}
\arguments{
  \item{df}{ snp.support data.frame }
  \item{tabname}{ name of SQLite table to be created and populated}
  \item{fn}{ filename of SQLite database store}
  \item{\dots}{ \dots, currently ignored}
}
\details{
This function makes an effort to ensure that chromosomes
are properly ordered in the resulting table.
The allele "assignments" from read.HapMap.data are translated
to IUPAC ambiguity codes to save space.
}
\value{
three data frames are invisibly returned:
the original input, the sorted data frame, and the
final data frame before serialization to SQLite
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
# would need a small support data frame and a
# temporary SQLite store
}
\keyword{ models }
