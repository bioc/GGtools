\name{make_racExSet}
\alias{make_racExSet}
\title{ create a racExSet from simpler constituents }
\description{create a racExSet from simpler constituents
}
\usage{
make_racExSet(exprs, racs, pd, mi, anno)
}
\arguments{
  \item{exprs}{ matrix of processed expression values, GxN,
 rownames are reporterNames}
  \item{racs}{ matrix of rare allele counts, SxN, rownames
are snp ids (rs numbers)}
  \item{pd}{ old phenoData instance}
  \item{mi}{ MIAME instance }
  \item{anno}{ string for annotation referencing}
}
%\details{
%}
%\value{
%}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(litsn)
data(litex)
data(pdd)
dem = make_racExSet( exprs=litex, racs=t(litsn), pd=pdd, mi=new("MIAME"),
  anno="hgfocus")
dem
}
\keyword{ models }
