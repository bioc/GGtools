
setClass("chunksize", contains="numeric")
chunksize = function(x) new("chunksize", as.numeric(x))


setMethod("gwSnpTests", c("formula", "smlSet", "snpdepth", "chunksize"),
 function(sym, sms, cnum, cs) {
# assumes a gene set is response of formula
  theCall = match.call(call=sys.call(2))
  gn = geneIds(gs <- eval(sym[[2]]))
  ng = length(gn)
   chunklabs = function (n, chunksize) 
   {
       bas = 1:n
       tool = ceiling(n/chunksize)
       as.numeric(cut(bas, tool))
   }
  gspl = split(gn, chunklabs(ng, cs))
  csets = lapply( gspl, function(x) gs[x] )
  savesym = sym
  out = list()
  for (i in 1:length(gspl)) {
      nsym = savesym
      nsym[[2]] = csets[[i]]
      out[[i]] = gwSnpTests(nsym, sms, cnum)
      gc()
      }
# this list has all the tests filtered already, so filterSnpTests is not 
# needed
  flattened = unlist(out, recursive=FALSE)
  names(flattened) = gn
  ans = new("filteredMultiGwSnpScreenResult", geneset=gs,
      call=theCall, flattened)
  names(ans@.Data) = gn
  ans
  })
   
