Strains2rac = function(strfile) {
#
# this function is intended to transform strain summary files from
#  http://www.well.ox.ac.uk/mouse/INBREDS into data structures
#  useful for GGtools Bioconductor package
# 
# note in the build34 strains file, there are some genotype calls H where all
# non-H strains are monomorphic -- can't determine the other base from these data
# they are set to NA
#
 tmp = read.table(strfile, h=TRUE, row.names=NULL)
 gt = t(t(tmp[,-c(1,2,3)]))
 alleles = function (x) 
  {
  # works for strains/gt format
    if (any(x == "HH")) 
        x = x[x != "HH"]
    sort(unique(unlist(strsplit(unique(x[x != "UU"]), ""))))
  }
 gt[] = paste(gt,gt, sep="")
 allelecodes = apply(gt,1,alleles)
 for (i in 1:nrow(gt)) if (length((AC <- allelecodes[[i]])>1)) 
      if (any(gt[i,]=="HH")) gt[i, gt[i,]=="HH"] = ifelse(length(AC)==2,paste(AC,collapse=""),NA)
 gt[gt == "UU"] = NA
 strains = names(tmp)[-c(1,2,3)]
 chr = as.character(tmp[,2])
 snp = as.character(tmp[,1])
 pos = as.numeric(as.character(tmp[,3]))
 rownames(gt) = snp
 colnames(gt) = strains
 rac = t(apply(gt, 1, countRare))
 rarebase =apply(gt,1, getRare)
 sla = function (x) 
   {
       if (length(x) == 1) 
           return(paste(x, x, sep = "/"))
       else if (length(x) == 2) {
           x = sort(x)
           return(paste(x[1], x[2], sep = "/"))
           }
       else return(NA)
   }
 snpAll = sapply(allelecodes, sla)
 list(gt=gt, chr=chr, pos=pos, rac=rac, rarebase=rarebase, SNPalleles=snpAll)
}

INBREDSworkflow = function(inbfile, emat, estrains, pd, mi, anno, fixup=NULL,
   fixchr=function(x)gsub("_random", "", x)) {
 srac = Strains2rac(inbfile)
 sn = estrains
 ssn = colnames(srac$gt)
 if (!is.null(fixup)) {
	ssn = fixup(ssn)
	colnames(srac$gt) = ssn
	}
 if (!all(sn %in% ssn)) stop(paste("there are some strain names in estrains \n[",
       paste(setdiff(sn,ssn),collapse=", "), "]\n that are not present in the inbfile columns\n",
"each expr sample strain must match to some inbfile column"))
 res = make_racExSet( emat, srac$gt[,estrains], srac$rarebase, srac$SNPalleles, pd, mi, anno )
 smw = wrapSNPmetaWh( rownames(srac$gt), fixchr(srac$chr), srac$pos )
 list(racExSet=res, snpMetaWh=smw)
}
 
