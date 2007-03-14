genoStrings = function(racExSet, rsnum) {
#
# this assumes that rarebase is named -- and it is not required to
# be.  so get an index
#
 if (length(rsnum)>1) stop("sorry, only handling a single snp id at present")
 wh = which(snpNames(racExSet) == rsnum)
 rac = snps(racExSet)[wh,]
 rar = rarebase(racExSet)[wh]
 gt = SNPalleles(racExSet)[wh]
 gts = as.character(strsplit(gt,"/")[[1]])
 com = gts[gts!=rar]
 opts = c( homrar=paste(rar,rar,sep="/"), het=as.character(gt), homcom=paste(com,com,sep="/"))
 crac = as.character(rac)
 crac[rac==0] = opts["homcom"]
 crac[rac==1] = opts["het"]
 crac[rac==2] = opts["homrar"]
 crac
}
 
