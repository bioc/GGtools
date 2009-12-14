\name{slimCTchr2SQLite}
\alias{slimCTchr2SQLite}
\alias{slimCTchr2flat}
\alias{slimCTchr2ff}
\alias{ffs2minp}
\alias{minp}
\alias{locreport}
\alias{wgtinfo2browser}
\alias{top4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transform matrices of chromosome-wide cis/trans eQTL statistics to
SQLite or flat text, and support extraction of minimum p per SNP,
with location information
}
\description{
transform matrices of chromosome-wide cis/trans eQTL statistics to
SQLite, flat text, or ff package managed
flat file format, and support extraction of minimum p per SNP
}
\usage{
slimCTchr2SQLite(path, chr = 1, dbname = paste("chr", chr, sep = ""))
slimCTchr2flat(path, chr = 1, dbname = paste("chr", chr, sep = ""))
slimCTchr2ff(path, chr = 1, dbname = paste("chr", chr, sep = ""))
minp(path, BUFSIZE=100000, breakat=Inf, freqdump=10000, filefun=gzfile)
ffs2minp(ffl, ...)
locreport(winfo, chr)
wgtinfo2browser(winfo, chr, start, end, trname="newtrack",
  bsession=NULL, vlim=c(0,16), ...)
top4(x, sms)
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
biological name of chromosome; for locreport or wgtinfo2browser, a numerical
chromosome identifier
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
\item{start}{starting base position for browser}
\item{end}{starting base position for browser}
\item{trname}{display name for new track on browser}
\item{bsession}{browserSession instance if you want to modify an existing session}
  bsession=NULL, ...)
\item{x}{list, result of splitting a locreport result by gene}
\item{sms}{an smlSet instance}
\item{vlim}{a 2-vector like the output of base::range() prescribing the
initial vertical range for the track to be browserd by wgtinfo2browser; units are
-log 10 p}
\item{ffl}{a list of ff references -- note that the ff package has to be
attached and you must load the .rda object corresponding to each ff entity}
\item{\dots}{for wgtinfo2browser, parameters passed to browserSession; for ffs2minp,
not currently used.}
}
\details{
uses dbWriteTable for SQLite driver built in function call,
or write.table (for 'flat' application)

note that wgtinfo2browser will winsorize 0 p-values to the smallest
non-zero values reported (on -log scale)
}
\value{
slim* NULL, just runs for side effects

minp -- a list of 3 elements, rsid (as numeric, no rs prefix),
ans (minimum p value), and hline (header record, enumerates
probe names)

wgtinfo2browser -- a ucsc browsersession instance
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
}
\note{
After limited experimentation it was decided that SQLite was not
effective for this application, while ff seems to be.

slimCTchr2ff takes a folder assumed to be populated with
lists of slimCisTrans results.  Each set of results is over
a common set of SNPs typically on multiple chromosomes, but on an
element of a partition of the set of genes.  slimCTchr2ff will
combine results over genes for all SNP on a single
chromosome and generate an ff array holding all
the information.  column 1 of the resulting array is the numeric
component of the dbSNP id for the SNP defining the row (it will
be NA if a non rs-number has been used), and the remaining
columns are snp-specific, gene-specific chisquared stats with 1df.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
overb = options()$verbose
example(slimCisTrans)
setwd(td)
#
# note that the chr=1 below indicates the first in sequence
# in this smlSet (hmceuB36.2021), which is chr20
#
\dontrun{
slimCTchr2flat(td, chr=1, dbname="chr20flat")
options(verbose=TRUE)
mp = minp("chr20flat.txt.gz")
sapply(mp, head, 10)
options(verbose=overb)
lr = locreport(mp, 20)
lr[1:10,]
# now visualize the best eQTL association on chr20
oo = order(lr[,"mlogp"], decreasing=TRUE)
#
# following throws margin error in CMD check for 2.10.0
#
if (interactive()) plot_EvG(probeId(lr[oo[1],"bestgene"]), rsid(
   lr[oo[1], "rsid"]), hmceuB36.2021[chrnum(20),])
lrc = na.omit(lr)
}
#
# now put the best scores for all snps on chr20 vs 10 genes on chr20
# on the browser
#
if (interactive()) {
 data(chr17wgtinfo)
 wgtinfo2browser(chr17wgtinfo, 17, 1, 20000000, "cis17")
}
#
# new interface uses ff structures instead of flat files
#
if (require(ff)) {
  slimCTchr2ff(td, chr=1, dbname="chr20flat")
  load("chr20flat.rda")
  chr20flatchisq
  mmf = ffs2minp( list(chr20flatchisq ) )
  sapply(mmf, length)
  sapply(mmf, class)
  sapply(mmf, head, 10)
  options(verbose=overb)
  lr = locreport(mmf, 20)
  lr[1:10,]
  oo = order(lr[,"mlogp"], decreasing=TRUE)
  if (interactive()) plot_EvG(probeId(lr[oo[1],"bestgene"]), rsid(
   lr[oo[1], "rsid"]), hmceuB36.2021[chrnum(20),])
  lrc = na.omit(lr)
#
# here we look at the best snps per gene, expression vs genotype,
# for four genes
#
  oo2 = order(lrc[,"mlogp"], decreasing=TRUE)
  slr = split(lrc[oo2,], lrc[oo2,"bestgene"])
  topsnps = sapply(slr, function(x)x[1,"mlogp"])
  olist = order(topsnps, decreasing=TRUE)
  top4(slr[olist[1]], hmceuB36.2021)
  }
# end examples
}
\keyword{ models }
