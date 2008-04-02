\name{plot_EvG-methods}
\docType{methods}
\alias{plot_EvG-methods}
\alias{plot_EvG}
\alias{plot_EvG,genesym,rsNum,smlSet-method}
\alias{getAlleles,smlSet,rsNum-method}
\alias{getAlleles}
\title{ plot expression vs. genotype }
\description{
plot expression vs. genotype ; getAlleles is supporting
infrastructure that may have additional uses
}
\section{Methods}{
\describe{
\item{gsym = "genesym", rsn = "rsNum", sms = "smlSet"}{ displays distributions
of gene expression values for gene
identified by \code{gsym} (which must be a symbol identified
in the [platform]SYMBOL AnnDbBimap for the 
smlSet annotation) stratified by SNP genotype  for SNP identified
by \code{rsn}.}
\item{getAlleles(sms="smlSet", rsn="rsNum")}{return character
tokens describing subject-specific genotypes}
}
}
\examples{
data(hmceuB36.2021)
plot_EvG( genesym("CPNE1"), rsNum("rs6060535"), hmceuB36.2021 )
getAlleles(hmceuB36.2021, rsNum("rs6060535"))[1:5]
}
\keyword{methods}
