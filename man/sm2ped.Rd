\name{sm2ped}
\alias{sm2ped}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a data.frame representing a PED file (MACH) from a SnpMatrix instance
}
\description{
create a data.frame representing a PED file (MACH) from a SnpMatrix instance
}
\usage{
sm2ped(sm, snpsupp, missing.code = "N", family, person, father, mother, sex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sm}{
\code{\link[snpStats]{SnpMatrix-class}} instance
}
  \item{snpsupp}{
data frame of 'supporting' metadata as returned by \code{\link[chopsticks]{read.HapMap.data}}
}
  \item{missing.code}{
token to use for missing genotype
}
  \item{family}{
vector of family identifiers as used in PED files for MACH
}
  \item{person}{
vector of person identifiers
}
  \item{father}{
each element of this vector specifies the father of the person identified on the
same line of the file
}
  \item{mother}{
as for \code{father}
}
  \item{sex}{
1 for male, 2 for female
}
}
\details{
basic purpose is to convert SnpMatrix instances to PED to allow genotype imputation
}
\value{
a data.frame instance
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(snpStats)
sm = chopsticks::read.HapMap.data(paste("file://",
 system.file("hapmap/smallc20CEU.txt.gz", package="GGtools"), sep=""))
supp = sm[[2]]
smat = sm[[1]]
reldata = read.table(system.file("hapmap/relationships_w_pops_051208.txt",
  package="GGtools"), header=TRUE)
rownames(reldata) = as.character(reldata[,2])
hids = rownames(smat)
hrel = reldata[hids,]
hrel[1:5,]
args(sm2ped)
ac = as.character
ss = sm2ped(smat, supp, fam=ac(hrel[,1]), person=ac(hrel[,2]), father=ac(hrel[,3]),
mother=ac(hrel[,4]), sex=hrel[,5])
ss[1:5,1:8]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
