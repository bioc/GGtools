GeneSet2LocInfo = function(gs) {
#
# returns a list of named vectors, each element has a gene set
# id as name and a vector of signed locations on chromosomes
# the names of the location elements are the chromosome names
#
 if (organism(gs) != "Homo sapiens") stop("only written for Homo sapiens at this time")
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
   toks = geneIds(respObj)
   lnames = lapply(geneLocList, names)
   lnl = lapply(lnames, nchar)
   for (i in 1:length(lnl))
     lnames[[i]] = lnames[[i]][ lnl[[i]] < 3 ]
   chroms = sapply(lnames, "[", 1)
   keepSnps = snpsNear(respObj, radius)
   ntests = length(chroms)
   conditions = list()
   outl = list()
   curfmla = fmla
   for (i in 1:ntests) {
     if (options()$verbose) cat(".")
     resp = wrapex( toks[i] )
     curfmla = formula(paste(resp , pred, sep="~"))
     cursm = smls[ chrnum(chroms[i]), ]
     smat = smList(cursm)[[chroms[i]]]
     onc = colnames(smat)
     nsnps = length(intersect(onc,keepSnps[[i]]) )
     if (nsnps == 0) {
         warning(paste("no snps on chip for given radius relative to gene; need to increase; executing full chromosome test for gene ", toks[i],"chr", chroms[i]))
         conditions[[i]] = list(gene=toks[i], chrom=chroms[i], cond="no SNP in radius")
         outl[[i]] = NA
       }
     else {
        smat = smat[, intersect(onc,keepSnps[[i]]) ]
        tmp = list(smat)
        names(tmp) = chroms[i]
        assign("smList", tmp, cursm@smlEnv)
        outl[[i]] = gwSnpTests(curfmla, cursm, chrnum(chroms[i]))
        }
     }
   return(new("multiCisTestResult", call=mc, conditions=conditions, outl))
  }
  else stop("only functioning for GeneSets based on probeIds")
}
  
