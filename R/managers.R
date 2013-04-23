

#
# an eqtlTestsManager can cover a collection of SNP on different
# chromosomes with a single set of genes
# fffile slot holds a single ff matrix (used to be list) where rows are SNP and columns are
#     genes
# call, sess, exdate geneanno slots are metadata
# shortfac is the scaling factor used to inflate chisq stats so short integer
#     representation has some precision on division by shortfac
# df is d.f. of chisq stat
#
# if em is an eqtlTestsManager instance then em[rsid, probeId] returns
#     a list of chisq statistics properly rescaled
#


chkeman = function(object){
# eqtlTestsManager validity test
 return(TRUE)
}

setMethod("shortfac", "eqtlTestsManager", function(x)
  x@shortfac)
setMethod("fffile", "eqtlTestsManager", function(x)
  x@fffile)
setMethod("exdate", "eqtlTestsManager", function(x)
  x@exdate)

setMethod("show", "eqtlTestsManager", function(object) {
 cat(class(object), " computed", exdate(object), "\n")
 cat("gene annotation:", object@geneanno, "\n")
 on.exit(close(fffile(object)))
 open(fffile(object))
 cat("some genes (out of ", length(colnames(fffile(object))),"): ", 
    paste(selectSome(colnames(fffile(object))),collapse=" "), "\n", sep="")
 cat("some snps (out of ", nrow(fffile(object)),  "): ", 
    paste(selectSome(rownames(fffile(object))),collapse=" "), "\n", sep="")
})

setMethod("[", c("eqtlTestsManager"), # , "rsid", "probeId"),
 function(x, i, j, ..., drop=FALSE) {
 on.exit(try(close(x@fffile)))  # new 13 may 2012 -- more consistent open/close for ff
 if (!is.open(x@fffile)) { # squelch ff warning by explicitly opening
   tmp = open(x@fffile)
   if (!is.open(x@fffile)) stop("failed to open ff archive")
   }
 sn = rownames(x@fffile)
 pn = colnames(x@fffile)
 if (!missing(i)) {
  if (!is(i, "numeric") && length(setdiff(i,sn))>0) {
    warning("scores for some nonexistent SNP were requested; these were dropped")
    i = intersect(i,sn)
    }
  }
 else i = ff::hi(1, length(sn))
 if (!missing(j)) {
  if (!is(j , "numeric") && length(setdiff(j,pn))>0) {
    warning("scores for some nonexistent probes were requested; these were dropped")
    j = intersect(j,pn)
    }
  }
 else j = ff::hi(1, length(pn))
 x@fffile[i,j, ..., drop=drop]/x@shortfac
})

probesManaged = function(x, ...) colnames(x@fffile)
snpsManaged = function(x, ...) rownames(x@fffile)

#setMethod("[", c("eqtlEstimatesManager"), # , "rsid", "probeId"),
# function(x, i, j, ..., drop=FALSE) {
# on.exit(try(close(x@fffile)))  # new 13 may 2012 -- more consistent open/close for ff
# if (!is.open(x@fffile)) { # squelch ff warning by explicitly opening
#   tmp = open(x@fffile)
#   if (!is.open(x@fffile)) stop("failed to open ff archive")
#   }
# sn = rownames(x@fffile)
# pn = colnames(x@fffile)
# if (!missing(i)) {
#  if (!is(i, "numeric") && length(setdiff(i,sn))>0) {
#    warning("scores for some nonexistent SNP were requested; these were dropped")
#    i = intersect(i,sn)
#    }
#  }
# else i = ff::hi(1, length(sn))
# if (!missing(j)) {
#  if (!is(j , "numeric") && length(setdiff(j,pn))>0) {
#    warning("scores for some nonexistent probes were requested; these were dropped")
#    j = intersect(j,pn)
#    }
#  }
# else j = ff::hi(1, length(pn))
# x@fffile[i,j, ..., drop=drop]/x@shortfac
#})
