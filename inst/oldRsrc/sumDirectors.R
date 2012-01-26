sumDirectors = function( dirl, prefix ) {
 if (missing(prefix)) stop("must supply prefix for ff files")
 mm = sapply(dirl, function(x) names(x@mgrs))
 if (!all(apply(mm,1,function(x)all(x[1]==x)))) stop("chrom names not equivalent between directors")
 chn = mm[,1]
 mg1 = dirl[[1]]@mgrs[[1]]
 alldf = sapply(dirl, function(x) x@mgrs[[1]]@df)
 if (!all(alldf == 1)) stop("assumes df = 1 for all managers") # superficially check first only in each
#
#Name:       fflist        call        sess      exdate    shortfac    geneanno
#Class:        list        call         ANY         ANY     numeric   character
#                              
#Name:           df summaryList
#Class:     numeric        list
  thecall = match.call()
  sess = sessionInfo()
  exdate = date()
  shortfac = mg1@shortfac
  geneanno = mg1@geneanno
  df = length(dirl)
  summaryList = list()
  nchr = length(dirl[[1]]@mgrs)
  dims = lapply(dirl[[1]]@mgrs, function(x)dim(x@fflist[[1]]))  # each manager assumed length 1 -- "cis"
  dimn = lapply(dirl[[1]]@mgrs, function(x)dimnames(x@fflist[[1]]))  # ditto
  ffl = list()
  mglist = list()
  for (i in 1:nchr) {
    ffl[[chn[i]]] = ff(0, vmode="short", filename=paste(prefix, "_", chn[i], ".ff", sep=""),
       dim = dims[[i]], dimnames=dimn[[i]])
    for (j in 1:length(dirl)) {
       ffl[[chn[i]]][,,add=TRUE] = dirl[[j]]@mgrs[[chn[i]]]@fflist[[1]][]
       }
    close(ffl[[chn[i]]])
    mglist[[chn[i]]] = new("eqtlTestsManager", fflist=ffl[i], 
           call=thecall, sess=sess, exdate=exdate, shortfac=shortfac,
	   geneanno=geneanno, df=df, summaryList=list())
    }
   new("multiCisDirector", mgrs=mglist)
}
 
