\name{appraise}
\alias{appraise}
\alias{bindgwava}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
appraisal for eQTL prediction models
}
\description{
appraisal for eQTL prediction models
}
\usage{
appraise(dtab, 
  discretize = TRUE, 
  reduceToSNP = TRUE, 
  prefix, 
  folder = paste0(prefix, "_APPROUT"), 
  discfmlas_in = GGtools:::.discfmlas.demo, 
txlist = list(
    distcats = function(x) {
     cut(x$mindist, c(-1, seq(0, 200001, 50000)))
     },
    fdrcats = function(x) {
     fdrfac = cut(x$fdr, c(-.01, .05, .1, .25, .5, 1.01))
     relevel(fdrfac, "(0.5,1.01]")
     },
    mafcats = function(x) {
     maffac = cut(x$MAF,c(-0.01,.05, .1, .25, .51))
     relevel(maffac, "(-0.01,0.05]")
     },
    caddcats = function(x){
     cut(x$PHRED, c(-.01, 5, seq(10, 30, 10 ), 60))
    }
    ),
    cutts = c(-0.01,seq(0.015,.12,.015),.15), 
    names2check= GGtools:::.standardNames, maxit=30,
     savePinfer=FALSE)
# bindgwava( gwavadt, eqdt )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dtab}{
data.table instance as created by transforming
cisRun to GRanges and then to data.table, and
then adding CADD PHRED scores if available.  If CADD
PHRED scores are not available, the default formulas
should not be used.
}
  \item{discretize}{
logical telling whether binning to factors defined in txlist
should be performed
}
  \item{reduceToSNP}{
logical telling whether ranges should be reduced to 
unique SNP and FDR recomputed
}
  \item{prefix}{
character atom used to prefix objects saved and folder for
result objects
}
  \item{folder}{
folder name suffix
}
  \item{discfmlas_in}{
named list of model formulae
}
  \item{txlist}{
named list of functions that are used to
bin certain quantitative features of SNP
}
  \item{cutts}{numeric vector of thresholds for tabulation and discrete calibration}
  \item{names2check}{if NULL, ignored; if a character vector, function
will fail unless \code{all(names2check \%in\% names(dtab)}
 }
  \item{maxit}{numeric passed to bigglm as control parameter for maximum
number of iterations to use in modeling gwas hit probabilities}
  \item{savePinfer}{logical specifying whether the
inferred probabilities of GWAS involvement are retained}
%  \item{gwavadt}{data.table including the GWAVA pathogenicity scores
% with columns 'tss', 'unmatched', 'region'}
%  \item{eqdt}{data.table from a cisRun output}
}
\details{
The appraise function wraps many tasks used to appraise eQTL
collections in terms of predictive capacity.  Details
will be provided.
}
\value{
A folder is opened and objects are written
representing the test set (data.table on SNPs on even chromosomes),
the coefficients of predictive models built on training set
(SNPs on odd chromosomes), coefficients of linear regressions
of binary test outcomes for calibrating the model on test data,
and ROC AUC measures.

bindgwava uses simple data.table operations with match to
add three columns to eqdt, gwava_tss, gwava_unmat, and gwava_regi
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
