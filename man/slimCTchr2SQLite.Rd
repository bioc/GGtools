\name{slimCTchr2SQLite}
\alias{slimCTchr2SQLite}
\alias{slimCTchr2flat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transform matrices of chromosome-wide cis/trans eQTL statistics to
SQLite or flat text
}
\description{
transform matrices of chromosome-wide cis/trans eQTL statistics to
SQLite or flat text
}
\usage{
slimCTchr2SQLite(path, chr = 1, dbname = paste("chr", chr, sep = ""))
slimCTchr2flat(path, chr = 1, dbname = paste("chr", chr, sep = ""))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
folder where slimCisTrans wrote the *+.rda files
}
  \item{chr}{
an index into the list of chromosome-specific results written
by slimCisTrans -- if you only ran slimCisTrans on a selected
chromosome, use 1 for this, then set the dbname to reflect the
biological name of chromosome
}

  \item{dbname}{
name of table and prefix of name of file to be written; in
case of 'flat' application, just interpreted as prefix
}
}
\details{
uses dbWriteTable for SQLite driver built in function call,
or write.table (for 'flat' application)
}
\value{
NULL, just runs for side effects
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
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
example(slimCisTrans)
slimCTchr2SQLite(td,chr=1, dbname="chr20")
setwd(td)
dr = dbDriver("SQLite")
con = dbConnect(dr, "chr20.db")
dbListTables(con)
dbGetQuery(con, "select * from chr20 limit 5")
# notice the fiddling with the column names in sqlite!
slimCTchr2flat(td,chr=1, dbname="chr20flat")
}
\keyword{ models }
