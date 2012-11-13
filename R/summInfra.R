setClass("sensiCisInput", representation(
   cisMgrFiles = "character",
   cisMgrProperties = "list",
   probeannopk = "character"))

setMethod("show", "sensiCisInput", function(object) {
   cat("eQTL sensitivity analysis inputs:\n")
   cat( length(object@cisMgrFiles), "runs to summarize.\n")
})

setClass("sensiCisOutput", representation(
   byGene = "GRanges", bySNP = "GRanges", tabAtFDRB="ANY",
    input = "sensiCisInput", thecall="call", fdrbound="numeric", sessionInfo="ANY"))
setMethod("show", "sensiCisOutput", function(object) {
  cat("eQTL sensitivity analysis results; genewise FDR bound was", object@fdrbound, ".\n")
  cat("There were", nrun <- length(object@input@cisMgrFiles), "runs.\n")
  cat("Selected results, ordered by yield:\n")
  if (nrun < 6) print(object@tabAtFDRB)
  else {
    hh = head(object@tabAtFDRB,3)
    ta = tail(object@tabAtFDRB,3)
    bb = rbind(as.matrix(hh), rep.int("...", ncol(hh)), as.matrix(ta))
    rownames(bb)[4] = "..."
    print(bb, quote=FALSE, right=TRUE)
    }
  cat("use @byGene and @bySNP to see specific tests; @tabAtFDRB for full summary.\n")
})

  

setGeneric("sensanal", function(object, fdrbound) standardGeneric("sensanal"))
setMethod("sensanal", c("sensiCisInput", "numeric"), function(object, fdrbound=0.05) {
   thecall = match.call()
   fns = object@cisMgrFiles
#   obl = lapply(fns, function(x) get(load(x)))
   summ = .summarize( props = object@cisMgrProperties,
      obs = fns, probeannopk = object@probeannopk, fdrbound=fdrbound )
   new("sensiCisOutput", byGene = summ$byGene, bySNP = summ$bySNP,
          tabAtFDRB = summ$tabAtFDRB,
          input=object, thecall=thecall, fdrbound=fdrbound, sessionInfo=sessionInfo() )})

.summarize = function(
   props = list( run1=c(rad=NA, excl=NA, maf=NA, nperm=NA, npc=NA) ),
   obs = c( somemcwbRefs=NA ),
   probeannopk = "illuminaHumanv1.db", fdrbound = 0.05 ) {
#
# infrastructure acquisition
#
 require(GGtools)
 require(IRanges)
 require(GenomicRanges)
 require(probeannopk, quietly=TRUE, character.only=TRUE)
#
# data frame converter from mcwBestCis
#
 adf = function (x) 
 {
     z = as.data.frame(fullreport(x))
     z$probe = rownames(z)
     rownames(z) = NULL
     z
 }
#
# deserialize all mcwBestCis instances named
#
 if (!is(obs, "character")) stop("obs should be a vector of strings naming mcwBestCis instances")
 bygene = lapply( obs, function(x) get(load(x)))
 allc = sapply(bygene, class)
 if (!all(allc==allc[1])) {
   print(table(allc))
   stop("deserialized instances not all mcwBestCis")
   }
#
# convert to data frames
#
 
 bygene = lapply( bygene, adf )

#
# propagate run properties to data frames
#

 for (i in 1:length(bygene) ) {
   bygene[[i]]$excl = props[[i]]["excl"]
   bygene[[i]]$maf = props[[i]]["maf"]
   bygene[[i]]$nperm = props[[i]]["nperm"]
   bygene[[i]]$npc = props[[i]]["npc"]
   }

#
# concatenate all runs
#

 full.bygene = do.call( rbind, bygene )

#
# obtain the best fdr over all runs, per gene
#

 full.bygene = full.bygene[ order(as.character(full.bygene$probe)), ]

 fbsf = split(full.bygene$fdr, full.bygene$probe )

 bestfdr = unlist(lapply( fbsf, function(x) rep(min(x), length(x)) ))

 full.bygene$bestfdr = bestfdr

#
# order report according to bestfdr over runs
#

 full.bygene2 = full.bygene[ order(full.bygene$bestfdr), ]

#
# add annotation
#
 sn = paste("chr", full.bygene2[,1], sep="")
 allp = full.bygene2[,"probe"]
 locsrc = get(paste(gsub(".db", "", probeannopk), "CHRLOC", sep=""))
 locendsrc = get(paste(gsub(".db", "", probeannopk), "CHRLOCEND", sep=""))
 symsrc = get(paste(gsub(".db", "", probeannopk), "SYMBOL", sep=""))
 st = sapply(mget(allp, locsrc), "[", 1)
 en = sapply(mget(allp, locendsrc), "[", 1)
 str = ifelse(st < 0, "-", "+")
#
# create GRanges by gene, start and end locations may depend on annotation versions
#
 drops = union(which(is.na(st)), which(is.na(en)))
 if (length(drops)>0)  {
    byGene_meta = GRanges(seqnames=sn[-drops], 
         IRanges(abs(st[-drops]),abs(en[-drops])), strand=str[-drops])
         values(byGene_meta) = DataFrame(full.bygene2[-drops,-c(1,2,3,4,5)])
    }
 else { byGene_meta = GRanges(seqnames=sn, 
         IRanges(abs(st),abs(en)), strand=str)
        values(byGene_meta) = DataFrame(full.bygene2[,-c(1,2,3,4,5)])
      }
 sym = sapply(mget(values(byGene_meta)$probe, symsrc), "[", 1)
 values(byGene_meta)$sym = sym
 byGene_cis = byGene_meta
 #save(byGene_cis, file="byGene_cis.rda")
#
# create GRanges by SNP
#
 bysnp = GRanges(seqnames=sn, IRanges(full.bygene2$snploc, width=1))#, values=full.bygene2[,-c(1,2,3,4,5,6,7)])
 values(bysnp) = full.bygene2[,-c(1,2,3,4,5,6)]
 values(bysnp)$snpid = as.character(values(bysnp)$snpid)
 bySNP_cis = bysnp
 spr = values(bySNP_cis)$probe
 sspr = sapply(mget(spr, symsrc), "[", 1)
 values(bySNP_cis)$sym = sspr
 bySNP_bestcis = bySNP_cis
 byGene_bestcis = byGene_cis
#
# create a flattened table as summary
#
 FTT <- with(full.bygene2, ftable(nperm, npc, maf, radiusUsed, excl, fdrltbnd=fdr <= fdrbound)) 
 DFTT = as.data.frame(FTT)
 tabAtFDRB = DFTT[DFTT[,7]>0 & DFTT[,6] == TRUE,]
 rownames(tabAtFDRB) = NULL
 colnames(tabAtFDRB) = c("nperm", "npc", "MAF", "radius", "exclRad", "fdrBound", "NgenesWeQTL")
 tabAtFDRB[,"fdrBound"] = fdrbound
 tabAtFDRB = tabAtFDRB[ order(tabAtFDRB[, "NgenesWeQTL"], decreasing=TRUE), ]
 list(bySNP=bySNP_bestcis, byGene=byGene_bestcis, tabAtFDRB=tabAtFDRB)
}
