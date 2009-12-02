\name{slimCTchr2SQLite}
\alias{slimCTchr2SQLite}
\alias{slimCTchr2flat}
\alias{minp}
\alias{locreport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transform matrices of chromosome-wide cis/trans eQTL statistics to
SQLite or flat text, and support extraction of minimum p per SNP,
with location information
}
\description{
transform matrices of chromosome-wide cis/trans eQTL statistics to
SQLite or flat text, and support extraction of minimum p per SNP
}
\usage{
slimCTchr2SQLite(path, chr = 1, dbname = paste("chr", chr, sep = ""))
slimCTchr2flat(path, chr = 1, dbname = paste("chr", chr, sep = ""))
minp(path, BUFSIZE=100000, breakat=Inf, freqdump=20000, filefun=gzfile)
locreport(winfo, chr)
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
\item{winfo}{output of minp (weighting info)}
\item{chr}{for locreport -- numerical chromosome identifier}
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
#
# note that the chr=1 below indicates the first in sequence
# in this smlSet (hmceuB36.2021), which is chr20
#
slimCTchr2flat(td,chr=1, dbname="chr20flat")
options(verbose=TRUE)
mp = minp("chr20flat.txt.gz")
sapply(mp, head, 10)
options(verbose=overb)
lr = locreport(mp, 20)
lr[1:10,]
# now visualize the best eQTL association on chr20
oo = order(lr[,"mlogp"], decreasing=TRUE)
plot_EvG(probeId(lr[oo[1],"bestgene"]), rsid(
   lr[oo[1], "rsid"]), hmceuB36.2021[chrnum(20),])
#
# now put the best scores for all snps on chr20 vs 10 genes on chr20
# on the browser
#
# it is pretty slow to get going, and it does not set ct20d to
# full view mode, but you can do that manually
#
if (interactive()) {
 library(rtracklayer)
 lrc = na.omit(lr)
 lrco = lrc[ order(lrc$reflocs),]
 rr = RangedData(IRanges(start=lrco[,2], end=lrco[,2]),
         score=lrco$mlogp, space="chr20")
 ee = browserSession("UCSC")
 tn = trackNames(ee)
 vp = ucscTrackModes(hide=setdiff(tn, c("ruler", "cytoband")), full="ct20d")
 ee[["ct20d"]] = rr
 ff = browserView(ee,
   GenomicRanges(start=1, end=62435964,chrom="chr20"), 
    track=c("ruler", "cytoBand"), full="ct20d")
 }
}
\keyword{ models }
