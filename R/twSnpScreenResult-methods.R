
setMethod("show", "twSnpScreenResult", function(object) {
   cat("twSnpScreenResult\n")
   cat("call: ")
   print(object@call)
   cat("Genes (selection):\n", selectSome(names(object)))
   cat("\nFirst fit object:\n")
   cat("---\n")
   show(object[[1]])
   cat("--- [there are", length(object)-1, "more]\n")
   })

