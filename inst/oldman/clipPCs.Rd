\name{clipPCs}
\alias{clipPCs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simple approach to removal of principal components from smlSet
}
\description{
simple approach to removal of principal components from smlSet
}
\usage{
clipPCs(smlSet, inds2drop, center=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
instance of \code{\link[GGBase]{smlSet-class}}
}
  \item{inds2drop}{
numeric vector of PCs to be eliminated
}
  \item{center}{logical passed to \code{\link[stats]{prcomp}}.}
}
\details{
uses SVD and zeroes out selected eigenvalues before reassembly
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an smlSet instance with transformed expression data
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (.Platform$OS.type != "windows") {  # our build system not removing folders...
if (!exists("hmceuB36.2021")) hmceuB36.2021 <- getSS("GGtools", c("20", "21"))
library(illuminaHumanv1.db)
g20 = get("20", revmap(illuminaHumanv1CHR))
g20 = intersect(g20, featureNames(hmceuB36.2021))[1:25]
hmc = clipPCs(hmceuB36.2021, 1:4)
hmc = hmc[probeId(g20),]
pcs = prcomp(t(exprs(hmceuB36.2021)))$x
hmr = hmceuB36.2021[ probeId(g20), ]
pData(hmr) = data.frame(pData(hmr), pcs[,1:4])
hmc
# files persist on certain windows systems?
if (file.exists("clipdem")) try(system("rmdir clipdem"))
if (file.exists("clipfmla")) try(system("rmdir clipfmla"))
if (file.exists("clipfmlaNOPC")) try(system("rmdir clipfmlaNOPC"))
f1 = eqtlTests(hmc[chrnum("20"),], ~male, targdir="clipdem")
f2 = eqtlTests(hmr[chrnum("20"),], ~male+PC1+PC2+PC3+PC4, targdir="clipfmla")
f3 = eqtlTests(hmr[chrnum("20"),], ~male, targdir="clipfmlaNOPC")
system("rm -rf clipdem")
system("rm -rf clipfmla")
system("rm -rf clipfmlaNOPC")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
