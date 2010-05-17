#
# in this file we introduce new formal class structures for
# collections of eqtl test runs
#

#
# an eqtlTestsManager can cover a collection of SNP on different
# chromosomes with a single set of genes
# fflist slot holds a list of ff matrices where rows are SNP and columns are
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
 allgn = lapply(fflist(object), colnames)
 n1 = allgn[[1]]
 chk = sapply(allgn[-1], function(x)all.equal(x,n1))
 if (!all(chk)) return("fflist colnames not common to all elements")
 if (is.null(names(fflist(object)))) return("fflist elements lack names")
 return(TRUE)
}

setClass("eqtlTestsManager",
 representation(fflist="list", call="call", sess="ANY",
	exdate="ANY", shortfac="numeric", geneanno="character", df="numeric"),
        validity=chkeman)


setGeneric("shortfac", function(x)standardGeneric("shortfac"))
setMethod("shortfac", "eqtlTestsManager", function(x)
  x@shortfac)
setGeneric("fflist", function(x)standardGeneric("fflist"))
setMethod("fflist", "eqtlTestsManager", function(x)
  x@fflist)
setGeneric("exdate", function(x)standardGeneric("exdate"))
setMethod("exdate", "eqtlTestsManager", function(x)
  x@exdate)

setMethod("show", "eqtlTestsManager", function(object) {
 cat("eqtlTools results manager, computed", exdate(object), "\n")
 cat("There are", length(fflist(object)), "chromosomes analyzed.\n")
})

snpIdList = function(x) lapply(fflist(x), rownames)
geneIds = function(x) colnames(fflist[[1]])
permuterm = function(l) {
 lens = sapply(l, length)
 eln = rep(names(l), lens)
 dat = unlist(l)
 names(eln) = dat
 eln
}
 
snpIdMap = function(ids, x) {
 silx = snpIdList(x)
 m1 = lapply(silx, function(y){ 
         ans = match(y, ids)
         names(ans) = y
         ans
         })
 m2 = permuterm(m1)
 names(m2) = unlist(silx)
 split(names(m2[ids]), m2[ids])
}

setMethod("[", c("eqtlTestsManager", "rsid", "probeId"),
 function(x, i, j, ..., drop=FALSE) {
#
# ultimately this may not be exposed, serving only for deep
# testing, because a director database may be required for every
# manager
#
 m1 = snpIdMap( as(i, "character"), x )
 ans = lapply(1:length(m1), function(i) fflist(x)[[names(m1)[i]]][ m1[[i]], 
    as(j, "character"), drop=FALSE]/shortfac(x))
 names(ans) = names(m1)
 ans
})
 
eqtlTests = function(smlSet, rhs=~1-1,
   runname="foo", targdir="foo", geneApply=lapply, chromApply=lapply,
   shortfac = 100, computeZ=FALSE, ... ) {
 theCall = match.call()
 sess = sessionInfo()
 fnhead = paste(targdir, "/", runname, "_", sep="")
 geneNames = featureNames(smlSet)
 chrNames = names(smList(smlSet))
 ngenes = length(geneNames)
 nchr = length(chrNames)
 system(paste("mkdir", targdir))
 cres = chromApply( chrNames, function(chr) {
   snpdata = smList(smlSet)[[chr]]
   #targff = paste( fnhead, "chr", chr, "_", "g", gene, ".ff" , sep="" )
   targff = paste( fnhead, "chr", chr, ".ff" , sep="" )
   snpnames = colnames(snpdata)
   nsnps = ncol(snpdata)
   store = ff( initdata=0, dim=c(nsnps, ngenes), dimnames=list(snpnames, geneNames), vmode="short",
                 filename = targff )
   geneApply( geneNames, function(gene) {
     ex = exprs(smlSet)[gene,]
     fmla = formula(paste("ex", paste(as.character(rhs),collapse=""), collapse=" "))
     numans = snp.rhs.tests(fmla, snp.data=snpdata, data=pData(smlSet), family="gaussian", ...)@chisq
     if (computeZ) {
       numans = sqrt(numans)
       signl = snp.rhs.estimates( fmla, snp.data=snpdata, data=pData(smlSet), family="gaussian", ... )
       bad = which(unlist(lapply(signl, is.null)))
       if (length(bad)>0) signl[bad] = list(beta=NA)
       ifelse(unlist(signl)>=0, 1, -1)
       numans = numans*signl
     }
     miss = is.na(numans)
     if (any(miss) & !computeZ) numans[which(miss)] = rchisq(length(which(miss)), 1)
     if (any(miss) & computeZ) numans[which(miss)] = rnorm(length(which(miss)))
     store[, gene, add=TRUE] = shortfac*numans
     NULL
     }) # end gene apply
  store
  })  # end chr apply
  names(cres) = chrNames
  exdate = date()
  new("eqtlTestsManager", fflist=cres, call=theCall, sess=sess, 
        exdate=exdate, shortfac=shortfac, geneanno=annotation(smlSet),
        df=1)
}

# director for group of managers

chkmgrs = function(object) {
   mcl = sapply(mgrs(object), class)
   chkc = sapply(mgrs(object), function(x) is(x, "eqtlTestsManager"))
   if (!all(chkc)) return("mgrs slot must only contain list of entities inheriting from eqtlTestsManager")
   annos = sapply(mgrs(object), function(x)x@geneanno)
   if (!all(annos==annos[1])) return("managers do not have identical gene annotation source")
   sids = lapply(mgrs(object), snpIdList)
   slchk = sapply(sids, function(x) all.equal(x, sids[[1]]))
   if (!all(sapply(slchk,isTRUE))) return("managers do not have identical SNP lists")
   return(TRUE)
}
   
setClass("cisTransDirector", representation(mgrs="list", indexdbname="character", 
   shortfac="numeric", snptabname="character", probetabname="character", probeanno="character"),
   validity=chkmgrs)

setGeneric("mgrs", function(x) standardGeneric("mgrs"))
setMethod("mgrs", "cisTransDirector", function(x) x@mgrs)

mkCisTransDirector = function(dl, indexdbname, snptabname, probetabname, probeanno, commonSNPs=TRUE) {
 if (length(grep("\\.sqlite$", indexdbname))==0) stop("indexdbname must have suffix .sqlite")
 cat("creating SQLite database...")
 cd = new("cisTransDirector", indexdbname=indexdbname, shortfac=shortfac(dl[[1]]), mgrs=dl,
     snptabname=snptabname, probetabname=probetabname, probeanno=probeanno)
 mkDirectorDb(cd, commonSNPs)
 cat("...done\n")
 cd
}

mkDirectorDb = function(cd, commonSNPs=TRUE) {
#
# this seems a painful approach because it build a data.frame in memory before writing
# you will probably need to use inserts, iterating through managers
#
 alls = lapply(mgrs(cd), function(x)lapply(fflist(x), rownames)) #
 allg = lapply(mgrs(cd), function(x) colnames(fflist(x)[[1]])) # assumes all chroms have common genes assayed
 sdf = NULL
 nsmgrs = length(cd)
 if (commonSNPs) {
   rsids = unlist(alls[[1]])
   cn = names(alls[[1]])
   chrs = rep(cn, sapply(alls[[1]], length))
   mgr = rep(1, length(rsids))
   sdf = data.frame(snpid=rsids, chr=chrs, mgr=as.integer(mgr),
       stringsAsFactors=FALSE)
 }
 else stop("only handling managers with common SNP fields")

 gdfs = lapply(1:length(allg), function(i) data.frame(probeid=allg[[i]],
        mgr=i))
 gdf = gdfs[[1]]
 if (length(gdfs)>1) for (i in 2:length(gdfs)) gdf = rbind(gdf, gdfs[[i]])

 library(RSQLite)
 drv = dbDriver("SQLite")
 con = dbConnect(drv, dbname=cd@indexdbname)
 dbWriteTable(con, cd@snptabname, sdf, row.names=FALSE)
 dbWriteTable(con, cd@probetabname, gdf, row.names=FALSE)
 dbDisconnect(con)
}


setMethod("show", "cisTransDirector", function(object) {
 cat("eqtlTools cisTransDirector instance.\n")
 cat("there are", length(mgrs(object)), "managers.\n")
 cat("First:\n")
 show(mgrs(object)[[1]])
})

