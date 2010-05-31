\name{mkCisTransDirector}
\alias{mkCisTransDirector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create an object that manages a collection of eqtlTestManagers
}
\description{
Create an object that manages a collection of eqtlTestManagers
}
\usage{
mkCisTransDirector(dl, indexdbname, snptabname, probetabname, probeanno, commonSNPs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dl}{
list of eqtlManager instances
}
  \item{indexdbname}{
scalar character used to distinguish the director
}
  \item{snptabname}{
name to be used for the index of snp to chromosomes
}
  \item{probetabname}{
name to be used for the index of probes to managers
}
  \item{probeanno}{
platform annotation package name, e.g., "illuminaHumanv1.db"
}
  \item{commonSNPs}{
logical indicating whether all managers cover the same collection of
SNPs 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
instance of cisTransDirector class
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu?
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# see example(eqtlTests)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
