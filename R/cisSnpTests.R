GeneSet2LocInfo = function(gs) {
#
# returns a list of named vectors, each element has a gene set
# id as name and a vector of signed locations on chromosomes
# the names of the location elements are the chromosome names
#
 if (annotate::organism(gs) == "") warning("organism missing from gene set object, assuming Homo sapiens")
 else if (annotate::organism(gs) != "Homo sapiens") stop("only functioning for gene sets satisfying organism(gs) %in% c('', 'Homo sapiens')")
 gst = geneIdType(gs)
 tags = geneIds(gs)
 if (gst@type == "Annotation") {
   anp = gst@annotation
   if (length(grep("\\.db", anp)) == 0) {
      prefix = anp
      anp = paste(anp, ".db", sep="")
      }
   else prefix = gsub("\\.db", "", anp)
   }
 else {
   anp = "org.Hs.eg.db"
   prefix = "org.Hs.eg"
   }
 require(anp, character.only=TRUE, quietly=TRUE)
 if (gst@type == "Symbol")
  tags = unlist(mget(tags, revmap(org.Hs.egSYMBOL)))
 mget(tags, get(paste(prefix, "CHRLOC", sep="")))
 }

 

cisSnpTests = function(fmla, smls, radius, ...) {
 mc = match.call()
 wrapg = function(x) paste("genesym(\"", x, "\")", sep="")
 wrapex = function(x) paste("probeId(\"", x, "\")", sep="")

 okRespClNames = c("GeneSet", "genesym", "probeId")
 respObj = eval(fmla[[2]])
 if (!(any(sapply(okRespClNames, function(c) is(respObj,c)))))
   stop("dependent variable must inherit from one of GeneSet, genesym or probeId")
 if (is(respObj, "GeneSet")) {
   flist = as.list(fmla)
   if (length(flist[[3]]) > 1) pred = paste(flist[[3]][-1], collapse="+")
     else if (length(flist[[3]]) == 1) pred = flist[[3]]

   geneLocList = GeneSet2LocInfo(respObj)
   isna = sapply(geneLocList, function(x)is.na(x[1]))
   isnull = sapply(geneLocList, function(x)length(x[1])==0)
   if (any(isna) | any(isnull)) {
       kill = union(which(isna), which(isnull))
       geneLocList = geneLocList[-kill]
       respObj = respObj[-kill]
       warning(paste(length(kill), "genes in response object dropped for lack of location info"))
       }
   toks = geneIds(respObj)
   lnames = lapply(geneLocList, names)
   lnl = lapply(lnames, nchar)
   for (i in 1:length(lnl))  # this is to get rid of chromosome annos with weird names
     lnames[[i]] = lnames[[i]][ lnl[[i]] < 3 ]
   chroms = sapply(lnames, "[", 1)
   keepSnps = snpsNear(respObj, radius)  # will have try-errors
   targInfo = lapply(keepSnps, function(x)try(attr(x, "target")))
   ntests = length(chroms)
   conditions = list()
   outl = list()
   curfmla = fmla
   for (i in 1:ntests) {
     if (options()$verbose) cat(".")
     if (inherits(keepSnps[[i]], "try-error") |   # try-error attribute gets lost
        length(grep("rror", keepSnps[[i]])>0)) {
        conditions[[i]] = list(gene=toks[i], chrom=
            chroms[i], cond="snpsNear fails for this gene, location probably ambiguous", targInfo=NA)
        outl[[i]] = NA
        next
        }
     resp = wrapex( toks[i] )
     curfmla = formula(paste(resp , pred, sep="~"))
     cursm = smls[ chrnum(chroms[i]), ,drop=FALSE]
     smat = smList(cursm)[[chroms[i]]]
     onc = colnames(smat)
     nsnps = length(intersect(onc,keepSnps[[i]]) )
     if (nsnps == 0) {
         warning(paste("no snps on chip for given radius relative to gene; need to increase; executing full chromosome test for gene ", toks[i],"chr", chroms[i]))
         conditions[[i]] = list(gene=toks[i], chrom=chroms[i], cond="no SNP in radius", targInfo=targInfo[[i]])
         outl[[i]] = NA
       }
     else {
        smat = smat[, intersect(onc,keepSnps[[i]]), drop=FALSE ]
        tmp = list(smat)
        names(tmp) = chroms[i]
        assign("smList", tmp, cursm@smlEnv)
        outl[[i]] = gwSnpTests(curfmla, cursm, chrnum(chroms[i]))
        conditions[[i]] = list(gene=toks[i], chrom=chroms[i], cond=NA, targInfo=targInfo[[i]])
        }
     }
   names(outl) = toks
   return(new("multiCisTestResult", call=mc, conditions=conditions, outl))
  }
  else stop("only functioning for GeneSets based on probeIds")
}
  
