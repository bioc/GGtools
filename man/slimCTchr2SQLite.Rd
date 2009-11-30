\name{slimCTchr2SQLite}
\alias{slimCTchr2SQLite}
\alias{slimCTchr2flat}
\alias{minp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transform matrices of chromosome-wide cis/trans eQTL statistics to
SQLite or flat text, and support extraction of minimum p per SNP
}
\description{
transform matrices of chromosome-wide cis/trans eQTL statistics to
SQLite or flat text, and support extraction of minimum p per SNP
}
\usage{
slimCTchr2SQLite(path, chr = 1, dbname = paste("chr", chr, sep = ""))
slimCTchr2flat(path, chr = 1, dbname = paste("chr", chr, sep = ""))
minp(path, BUFSIZE=100000, breakat=Inf, freqdump=20000, filefun=gzfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
for minp, the path of the file holding flat output, otherwise
folder where slimCisTrans wrote the *+.rda files.  for minp
the file will typically be compressed; if not, set filefun=file (literally)
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
\item{BUFSIZE}{ numeric specification of a buffer size to
avoid incremental growth of arrays }
\item{breakat}{ for testing, line number at which to stop
prematurely}
\item{freqdump}{ rate at which a record number will be printed,
if options()$verbose == TRUE }
\item{filefun}{ R function to open connection to flat file with
data; for recent runs of slimCTchr2flat, gzfile will work because
the output file is gzcompressed.}
}
\details{
uses dbWriteTable for SQLite driver built in function call,
or write.table (for 'flat' application)
}
\value{
slim* NULL, just runs for side effects

minp -- a list of 3 elements, rsid (as numeric, no rs prefix),
ans (minimum p value), and hline (header record, enumerates
probe names)
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
overb = options()$verbose
example(slimCisTrans)
if (FALSE) {
slimCTchr2SQLite(td,chr=1, dbname="chr20")
setwd(td)
dr = dbDriver("SQLite")
con = dbConnect(dr, "chr20.db")
dbListTables(con)
dbGetQuery(con, "select * from chr20 limit 5")
# notice the fiddling with the column names in sqlite!
}
setwd(td)
slimCTchr2flat(td,chr=1, dbname="chr20flat")
options(verbose=TRUE)
mp = minp("chr20flat.txt.gz")
sapply(mp, head, 30)
options(verbose=overb)
}
\keyword{ models }
