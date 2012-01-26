
topSingSnps = function(x, n=250, df=1) {
  pp = p.value(x) #, df)
  opp = order(pp)
  x[opp[1:n]]
}

filterGWS = function(x, ...) {
  gc()
  x@.Data = lapply(x@.Data, topSingSnps, ...)
  new("filteredGwSnpScreenResult", x)
}


