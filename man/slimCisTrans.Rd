\name{slimCisTrans}
\alias{slimCisTrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lightweight output for cis-trans eQTL search
}
\description{
lightweight output for cis-trans eQTL search
}
\usage{
slimCisTrans(smlSet, genes2do = 1:45, ncores = 15, targdir = "/mnt/data/stvjc/GWAS", fmla = gs ~ male)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
instance of smlSet class
}
  \item{genes2do}{
numerical indices of genes to be analyzed
}
  \item{ncores}{
number of cores to be used in mclapply
}
  \item{targdir}{
folder where results are to be written
}
  \item{fmla}{
general modeling formula as used in
gwSnpTests, dependent variable must be gs
}
}
\details{
This function uses multicore library to distribute work.
You may see problems if you have some annotation packages
loaded on head that are also used by slaves.
}
\value{
Returns list of NULL.  Main effect is creation of
files of lists.  See example.
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%
%}
\examples{
require(multicore)
data(hmceuB36.2021)
td = tempdir()
m1 = slimCisTrans(fmla=gs~male, smlSet=hmceuB36.2021,
   genes2do=1:10, ncores=5, targdir=td)
dir(td)
od = getwd()
setwd(td)
load(dir(td)[1])
names(tobj <- get(gsub(".rda", "", dir(td)[1])))
tobj$chisq[[1]][1:3,]
setwd(od)
}
\keyword{ models }
