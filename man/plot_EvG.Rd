\name{plot_EvG}
\alias{plot_EvG}
\title{plot expression vs genotype  }
\description{plot expression vs genotype
}
\usage{
plot_EvG(reset, gene, snpid, anno="hgfocus", jitfac=0.5, ...)
}
\arguments{
  \item{reset}{ instance of \code{\link{racExSet-class}}}
  \item{gene}{ instance of \code{\link{genesym-class}}}
  \item{snpid}{ character string with RS number for SNP}
  \item{anno}{ character annotation stub, such as
\code{hgfocus}}
  \item{jitfac}{numeric argument passed to \code{jitter} for
displacing abscissae of points}
  \item{\dots}{additional parameters to be passed to \code{plot};
at present, \code{ylab, xlab, pch, axes, cex}, and \code{cex.lab} parameters
are hardcoded, so it is an error to pass these in this call.}
% plot(jitter(x,.5), y, ylab = paste("log", gene, "expression"),
%        xlab = paste("minor allele count,", snpid), pch = 20,
%        axes = FALSE, cex = 1.9, cex.lab=1.5, ...)
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
