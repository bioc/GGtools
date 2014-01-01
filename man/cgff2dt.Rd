\name{cgff2dt}
\alias{cgff2dt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
translate the GFF3 from a ciseqByCluster/processgff output into
a serialized data.table instance, compute genome-wide plug-in FDR,
and update the GFF3 with this FDR
}
\description{
translate the GFF3 from a ciseqByCluster/processgff output into
a serialized data.table instance, compute genome-wide plug-in FDR,
and update the GFF3 with this FDR
}
\usage{
cgff2dt(gff3, tiling)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gff3}{
character string naming a tabix-indexed, bgzipped output of
\code{\link{gffprocess}}
}
  \item{tiling}{
output of \code{\link{simpleTiling}}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
assumes unix utilites zcat, paste and bgzip are available
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (gff3, tiling) 
{
    require(foreach)
    require(Rsamtools)
    stopifnot(is(tiling, "GRanges"))
    basen = gsub(".gff3.gz", "", gff3)
    th = headerTabix(gff3)
    orderedChr = th$seqnames
    lgr = foreach(i = 1:length(tiling)) \%dopar\% {
        gc()
        cat(i)
        lk = try(import.gff3(gff3, which = tiling[i]))
        if (inherits(lk, "try-error")) 
            lk = NULL
        if (!is.null(lk)) 
            lk = as.data.table(as.data.frame(lk))
        lk
    }
    bad = sapply(lgr, is.null)
    if (any(bad)) 
        lgr = lgr[-which(bad)]
    ans = do.call(rbind, lgr)
    ans$snplocs = as.numeric(ans$snplocs)
    ans$ests = as.numeric(ans$ests)
    ans$se = as.numeric(ans$se)
    ans$oldfdr = as.numeric(ans$fdr)
    ans$MAF = as.numeric(ans$MAF)
    ans$dist.mid = as.numeric(ans$dist.mid)
    nperm = length(grep("permS", names(ans)))
    pnames = paste("permScore_", 1:nperm, sep = "")
    for (i in 1:nperm) ans[[pnames[i]]] = as.numeric(ans[[pnames[i]]])
    ans$mindist = as.numeric(ans$mindist)
    print(date())
    ans$fdr = pifdr(ans$score, c(ans$permScore_1, ans$permScore_2, 
        ans$permScore_3))
    print(date())
    obn = paste0(basen, "_dt")
    assign(obn, ans)
    save(list = obn, file = paste0(obn, ".rda"))
    gwfdr = ans$fdr
    gwch = ans$seqnames
    gwsnp = ans$snp
    gwprobeid = ans$probeid
    sgwfdr = split(gwfdr, gwch)
    sgwfdr = unlist(sgwfdr[orderedChr])
    nn = tempfile()
    fdrt = file(nn, "w")
    writeLines(c(c(" ", " ", " "), paste("; gwfdr=", round(sgwfdr, 
        6), sep = "")), con = fdrt)
    close(fdrt)
    chkb = system(paste0("zcat ", gff3, " | paste -d '' - ", 
        nn, " | bgzip > ", gsub(".gff3", "wmlf.gff3", gff3)))
    if (!(chkb == 0)) 
        warning(paste("final system zcat-paste-bgzip returned non-zero: chkb=", 
            chkb))
    invisible(chkb)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
