\name{best.trans.eQTLs}
\alias{best.trans.eQTLs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
collect strongest trans SNP-gene associations in a buffer of size K genes per SNP
}
\description{
collect strongest trans SNP-gene associations in a buffer of size K genes per SNP
}
\usage{
best.trans.eQTLs(smpack, rhs, genechrnum, snpchrnum, K = 20, 
   targdirpref = "tsco", batchsize = 200, radius = 2e+06, genequeryprefix = "", 
   snploadprefix = "chr", snplocprefix = "chr", geneannopk, snpannopk, 
   exFilter = function(x) x, smFilter = function(x) x, geneApply = lapply)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smpack}{
character string naming a package from which \code{\link[GGBase]{smlSet-class}} instances
can be generated using \code{\link[GGBase]{getSS}}
}
  \item{rhs}{
passed to \code{\link[snpStats]{snp.rhs.tests}} for covariate or stratification adjustments;
for permutation analysis, covariates should be handled via \code{\link{regressOut}}
}
  \item{genechrnum}{
character vector of chromosome identifiers for genes, typically \code{as.character(1:22)} for
somatic genes in human studies
}
  \item{snpchrnum}{
specific chromosome identifier for all SNP to be analyzed
}
  \item{K}{
the size of the buffer: scores will be recorded for the most strongly associated K genes for each SNP 
}
  \item{targdirpref}{
character string where buffer data will be held in ff archives
}
  \item{batchsize}{
passed to \code{\link[ff]{ffrowapply}} as scores are filtered from comprehensive testing to fill
the buffer
}
  \item{radius}{
numeric: for same-chromosome tests, tests will not be performed for SNP-gene combinations with
base-pair proximity smaller than radius
}
  \item{genequeryprefix}{
string: used when the numeric chromosome identifier requires a prefix like `chr' for annotation
query resolution on gene location
}
  \item{snploadprefix}{
string: used when the package identified in \code{smpack} requires a prefix to the \code{snpchrnum}
token for \code{getSS} retrieval of \code{smlSet} instance
}
  \item{snplocprefix}{
string: used when the numeric chromosome identifier requires a prefix like `chr' for annotation
query resolution on SNP location
}
  \item{geneannopk}{
package to be used for CHRLOC and CHRLOCEND queries for genes
}
  \item{snpannopk}{
package to be used to resolve getSNPlocs calls
}
  \item{exFilter}{
function returning an smlSet instance, operating on expression component prior to smFilter application and eQTL testing
}
  \item{smFilter}{
function returning an smlSet instance, operating on the full smlSet
}
  \item{geneApply}{
lapply-like function, typically mclapply or the like
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
instance of \code{\link{transManager-class}}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (.Platform$OS.type != "windows") {  # ff overwrites failing 5.IX.12
  nsFilter2 = function(sms, var.cutoff=.5) {
   alliq = apply(exprs(sms),1,IQR)
   qs = quantile(alliq,var.cutoff, na.rm=TRUE)
   sms[ which(alliq > qs), ]
  }
  thefilt = function(x) GTFfilter( nsFilter2 (clipPCs(x, 1:10), var.cutoff=.95 ), lower=.05 )
  tfile = tempfile()
  tfold = dir.create(tfile)
  t1 = best.trans.eQTLs( "GGdata", ~1, as.character(20:22), "22",
          geneannopk="illuminaHumanv1.db", snpannopk= "SNPlocs.Hsapiens.dbSNP.20111119",
          smFilter=thefilt, snploadprefix="", snplocprefix="ch", targdirpref=tfile)
  tt1 = transTab(t1)
  tt1o = tt1[ order(tt1[,"sumchisq"], decreasing=TRUE), ][1:10,]
  tt1o
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
