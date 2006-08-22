\name{make_racExSet}
\alias{make_racExSet}
\title{ create a racExSet from simpler constituents }
\description{create a racExSet from simpler constituents
}
\usage{
make_racExSet(exprs, racs, rarebase, SNPalleles, pd, mi, anno)
}
\arguments{
  \item{exprs}{ matrix of processed expression values, GxN,
 rownames are reporterNames}
  \item{racs}{ matrix of rare allele counts, SxN, rownames
are snp ids (rs numbers)}
  \item{rarebase}{ character S-vector of rare base }
  \item{SNPalleles}{character S-vector of strings of form "A/B" indicating
biallelic polymorphism options}
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
data(rb20M)
data(Sa20M)
dem = make_racExSet( exprs=litex, racs=t(litsn), pd=pdd, mi=new("MIAME"),
  rarebase=rb20M, SNPalleles=Sa20M, anno="hgfocus")
dem
}
\keyword{ models }
