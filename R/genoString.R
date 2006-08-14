genoStrings = function(racExSet, rsnum) {
 rac = snps(racExSet)[rsnum,]
 rar = rarebase(racExSet)[rsnum]
 gt = SNPalleles(racExSet)[rsnum]
 gts = as.character(strsplit(gt,"/")[[1]])
 com = gts[gts!=rar]
 opts = c( homrar=paste(rar,rar,sep="/"), het=as.character(gt), homcom=paste(com,com,sep="/"))
 crac = as.character(rac)
 crac[rac==0] = opts["homcom"]
 crac[rac==1] = opts["het"]
 crac[rac==2] = opts["homrar"]
 crac
}
 
