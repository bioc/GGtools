#setClass("extSmlSet", representation=representation(
#   packname="character", chrnames="character"))
#
#setMethod("show", "extSmlSet", function(object) {
#  cat("externalized smlSet: package", object@packname, "\n", sep=" ")
#  cat("chromosomes:", selectSome(object@chrnames), "\n", sep=" ")
#  cat("use getSS('", object@packname, "', [chrnamevec]) to retrieve.\n", sep="")
#})

#lit = getSS( "ceuhm2", "chr22" )
#fn =featureNames(lit)
#rm(lit)
#gc()
#    
#library(illuminaHumanv1.db)
#litmap = lapply(c("1", "2", "3"), function(x) get(x, revmap(illuminaHumanv1CHR)))
#litmap = lapply(litmap, function(x) intersect(x, fn)[1:200])
#names(litmap) = paste("chr", c("1", "2", "3"), sep="")
#
##set.seed(1234)
##system("rm -rf foo")
##unix.time(md <- makeDiagDirector( "ceuhm2", litmap ))
#library(multicore)
#options(cores=12)
#set.seed(1234)
#system("rm -rf foo2")
#unix.time(md2 <- makeDiagDirector( "ceuhm2", litmap , geneApply=mclapply, targdir="foo2"))

externalize = function(smlSet, packname, author="Replace Me <auth@a.b.com>",
  maintainer="Replace Me <repl@a.b.com>") {
# creates folder structure for package
# saves expression data as ex in data folder 
 system(paste("mkdir", packname))
 system(paste("mkdir ", packname, "/inst", sep=""))
 system(paste("mkdir ", packname, "/R", sep=""))
 system(paste("mkdir ", packname, "/inst/parts", sep=""))
 cn = names(smList(smlSet))
 partfol = paste(packname, "inst/parts/", sep="/")
 datfol = paste(packname, "data/", sep="/")
 for (i in cn) { assign(i, smList(smlSet)[[i]]); save(list=i, file=paste(partfol, i, ".rda", sep="")) }
 system(paste("mkdir ", packname, "/data", sep=""))
 ex = as(smlSet, "ExpressionSet")
 save(ex, file=paste(datfol, "eset.rda", sep=""))
 dd = readLines(system.file("extpacksupp/DESCRIPTION.proto", package="GGtools"))
 zz = readLines(system.file("extpacksupp/zzz.R", package="GGtools"))
 dd = gsub("@MAINTAINER@", maintainer, dd)
 dd = gsub("@AUTHOR@", author, dd)
 dd = gsub("@PKGNAME@", packname, dd)
 writeLines(dd, paste(packname, "/DESCRIPTION", sep=""))
 writeLines(zz, paste(packname, "/R/zzz.R", sep=""))
 writeLines("", paste(packname, "/NAMESPACE", sep=""))
 cat(paste("now install", packname, "\n"))
 NULL
}

getSS = function( packname, chrs ) {
 require(packname, character.only=TRUE)
 ex = get(load(dir(system.file(package=packname, "data"), full=TRUE)))
 partsfol = system.file("parts", package=packname)
 sml = lapply(chrs, function(x) get(load(paste(partsfol, "/", x, ".rda", sep=""))))
 names(sml) = chrs
 make_smlSet( ex, sml )
}

#hmceuB36 = new("extSmlSet", packname="ceuhm2", chrnames=
#  paste("chr", c(1:22,"X", "Y"), sep=""))

.makeDiagDirector = function(packname, genemap, rhs=~1, ...) {
#
# NB -- this worked a few times on rex but then died mysteriously repeatedly
# with camp data -- now trying to refrain from writing to a single ff from two cores...
#
# packname is a package of components of an smlSet made by externalize() and by hand
# genemap is a list with elements corresponding to chromosomes, each element is a vector of probe ids
# we construct a multiCisDirector with all same-chromosome tests
# ... is passed to eqtlTests
#
 require(packname, character.only=TRUE)
 cnames = gsub(".rda", "", dir(system.file("parts", package=packname)))
 gmnames = names(genemap)
 mgrs = list()
 if (!all(gmnames %in% cnames)) stop("some chr in gene map is not represented in chroms of package")
 for (i in 1:length(genemap)) {
    tmp = getSS( packname, gmnames[i] )
    tmp = tmp[ probeId( genemap[[i]] ), ]
    mgrs[[ gmnames[i] ]] = eqtlTests( tmp, rhs, ... )
 }
 new("multiCisDirector", mgrs = mgrs )
}

makeDiagDirector = function(packname, genemap, rhs=~1, geneApply=lapply,
    mapApply=lapply, ...) {
#
# packname is a package of components of an smlSet made by externalize() and by hand
# genemap is a list with elements corresponding to chromosomes, each element is a vector of probe ids
# we construct a multiCisDirector with all same-chromosome tests
# ... is passed to eqtlTests
#
 require(packname, character.only=TRUE)
 cnames = gsub(".rda", "", dir(system.file("parts", package=packname)))
 gmnames = names(genemap)
# mgrs = list()
 if (!all(gmnames %in% cnames)) stop("some chr in gene map is not represented in chroms of package")
# for (i in 1:length(genemap)) {
 mgrs = mapApply( 1:length(genemap), function(i) {
    tmp = getSS( packname, gmnames[i] )
    tmp = tmp[ probeId( genemap[[i]] ), ]
    eqtlTests( tmp, rhs, geneApply=geneApply, ... )
 })
 names(mgrs) = names(genemap)
 new("multiCisDirector", mgrs = mgrs )
}

makeMultiDiagDirector = function(packnames, genemap, rhslist=list(~1), mapapply=lapply, geneApply=lapply, ...) {
#
# it is assumed that all smlSets externalized in packnames are conformant -- same featureNames,
#    same snp sets.  failure of this assumption could lead to disaster
#
# packnames is a vector of names of
#     packages of components of smlSets made by externalize() and by hand
#
# genemap is a list with elements corresponding to chromosomes, 
#     each element is a vector of probe ids
#
# mapapply tells how to iterate over the map, each step computes tests
#     for all probes in each map element (chromosome), using geneApply
#     mclapply should work here
#
# geneApply tells how to iterate over genes to add information to the ff
#     constructed for each map element.  probably lapply is safest here
#
# we construct a multiCisDirector with all same-chromosome tests
# ... is passed to eqtlTests
#
 lapply(packnames, function(z) require(z, character.only=TRUE))
 cnames = gsub(".rda", "", dir(system.file("parts", package=packnames[1])))
 gmnames = names(genemap)
# mgrs = list()
 if (!all(gmnames %in% cnames)) stop("some chr in gene map is not represented in chroms of package")
 mgrs = mapapply( 1:length(genemap), function(i) {
    tmp = getSS( packnames[1], gmnames[i] )
    tmp = tmp[ probeId( genemap[[i]] ), ]
    thismgr = eqtlTests( tmp, rhslist[[1]], ... )
    rm(tmp)
    gc()
    for (j in 2:length(packnames)) {
      tmp = getSS( packnames[j], gmnames[i] )
      tmp = tmp[ probeId( genemap[[i]] ), ]
#
# following will assume ffind=1 ... default
#
      increment1( thismgr@fflist[[1]], 
              eqtlTestsNofile( tmp, rhslist[[j]], ... )@fflist[[1]] )
      rm(tmp)
      gc()
      }
  })
 new("multiCisDirector", mgrs = mgrs )
}

increment1 = function( ff1, ff2, nchunk=20 ) {
 nc = ncol(ff1)
 if (nc != ncol(ff2)) stop("incompatible args")
 nr = nrow(ff1)
 if (nr != nrow(ff2)) stop("incompatible args")
 inds = chunk( 1, nc, length.out=nchunk )
 for (i in 1:length(inds) )
   ff1[, inds[[i]], add=TRUE] = ff2[, inds[[i]] ]
 NULL
}
