\name{ciseqByCluster}
\alias{ciseqByCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
end-to-end cluster-based cis-eQTL search
}
\description{
end-to-end cluster-based cis-eQTL search
}
\usage{
ciseqByCluster(cl, pack = "yri1kgv", outprefix = "yrirun", chromsToRun = 1:22, 
  targetfolder = "/freshdata/YRI_3", radius = 100000L, nperm = 3L, ncoresPerNode = 8, 
  numPCtoFilter = 10, lowerMAF = 0.02, geneannopk = "lumiHumanAll.db", 
  snpannopk = "SNPlocs.Hsapiens.dbSNP.20120608", smchrpref = "chr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{
instance of S3 cluster class from parallel package
}
  \item{pack}{
character string naming package to which \code{\link{getSS}} can be
applied to generate \code{\link[GGBase]{smlSet-class}} instances
}
  \item{outprefix}{
character string used to prefix names of output GFF3 files
}
  \item{chromsToRun}{
numeric tags of chromosomes to be analyzed
}
  \item{targetfolder}{
character string naming folder where GFF3 will be deposited
}
  \item{radius}{
extent of search around gene model in bp
}
  \item{nperm}{
number of permutations for plug-in FDR computation (usually a small integer)
}
  \item{ncoresPerNode}{
number of cores for multicore testing: chromosomes map to nodes, genes map to cores
}
  \item{numPCtoFilter}{
number of PCs to be removed through \code{\link{clipPCs}}
}
  \item{lowerMAF}{
lower bound on MAF of SNP to be included for testing
}
  \item{geneannopk}{
character string naming Bioconductor package with annotation for expression probe identifiers
}
  \item{snpannopk}{
character string naming Bioconductor package with annotation for SNP locations
}
  \item{smchrpref}{
character prefix converting chromsToRun elements to basenames of rda files  harboring SnpMatrix instances
}
}
\details{
purpose is to maximize throughput of cis-eQTL testing in a two-level concurrent computing environment,
where a cluster as defined in package parallel has nodes to which half-chromosomes 
will be dispatched; each
node is assumed to be multicore and genes are mapped to cores during the iteration process.
}
\value{
a set of GFF3 files encoding all cis associations with location and various metadata
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gffprocess}}, \code{\link{cgff2dt}}
}
\examples{
#none yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
