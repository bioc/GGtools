\name{plot_EvG}
\alias{plot_EvG}
\title{plot expression vs genotype  }
\description{plot expression vs genotype
}
\usage{
plot_EvG(reset, gene, snpid, anno="hgfocus")
}
\arguments{
  \item{reset}{ instance of \code{\link{racExSet-class}}}
  \item{gene}{ instance of \code{\link{genesym-class}}}
  \item{snpid}{ character string with RS number for SNP}
  \item{anno}{ character annotation stub, such as
\code{hgfocus}}
}
\details{
renders a simple scatter plot of
expression vs rare allele count
}
\value{
}
\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{  }


\seealso{  }

\examples{
data(chr20GGdem)
plot_EvG( chr20GGdem, genesym("CPNE1"), "rs6060535",
  "hgfocus")
}
\keyword{ models }
