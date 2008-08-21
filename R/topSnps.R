
setGeneric("topSnps", function(x, ...) standardGeneric("topSnps"))
setMethod("topSnps", "cwSnpScreenResult", function(x, n=10, which="p.1df") {
   DF = 1
   if (which == "p.2df") DF = 2
   else if (which  != "p.1df") warning("chisq df assumed to be 1, 'which' not recognized")
   pp = p.value(x@.Data[[1]], DF)
   spp = pp[ order(pp, decreasing=FALSE) ]
   df = data.frame(p.val=spp)
   rownames(df) = names(spp)
   df[1:n,,drop=FALSE]
})

setMethod("topSnps", "gwSnpScreenResult", function(x, n=10, which="p.1df") {
  ts.df = function (w, n = 10, which = "p.1df") {
   DF = 1
   if (which == "p.2df") DF = 2
   else if (which  != "p.1df") warning("chisq df assumed to be 1, 'which' not recognized")
   pp = p.value(w, DF)
   spp = pp[ order(pp, decreasing=FALSE) ]
   df = data.frame(p.val=spp)
   rownames(df) = names(spp)
   df[1:n,,drop=FALSE]
   }
  lapply(x, ts.df, n=n, which=which)
})

