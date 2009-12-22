\name{makeCommonSNPs}
\alias{makeCommonSNPs}
\alias{checkCommonSNPs}
\title{ confine the SNPs (in multiple chromosomes) in all elements of
a list of smlSets to the largest shared subset per chromosome;
test for satisfaction of this condition }
\description{
confine the SNPs (in multiple chromosomes) in all elements of
a list of smlSets to the largest shared subset per chromosome;
test for satisfaction of this condition
}
\usage{
makeCommonSNPs(listOfSms)
checkCommonSNPs(listOfSms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listOfSms}{
an R list with each element consisting of a
\code{\link{smlSet-class}}
}
}
\details{
intersection of set of rsids per chromosome is computed over all elements
}
\value{
list of smlSet instances sharing all SNP on all chromosomes
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(smlSet.example)
tmp = smList(smlSet.example)[[1]]
tmp = tmp[,-c(20:40)]
newe = new.env()
assign("smList", list(`21`=tmp), newe)
ex2 = smlSet.example
ex2@smlEnv = newe
try(checkCommonSNPs(list(smlSet.example,ex2)))
list2 = makeCommonSNPs( list(smlSet.example, ex2) )
checkCommonSNPs(list2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
